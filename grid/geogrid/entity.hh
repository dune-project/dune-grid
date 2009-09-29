// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/referenceelements.hh>

#include <dune/grid/geogrid/capabilities.hh>
#include <dune/grid/geogrid/cornerstorage.hh>

namespace Dune
{


  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class Grid,
      bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v) >
  class GeometryGridEntityImpl;

  template< int codim, int dim, class Grid >
  class GeometryGridEntity;



  // GeometryGridEntityImpl (real)
  // -----------------------------

  template< int codim, class Grid >
  class GeometryGridEntityImpl< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

  public:
    typedef typename Traits :: ctype ctype;

    enum { codimension = codim };
    enum { dimension = Traits :: dimension };
    enum { mydimension = dimension - codimension };
    enum { dimensionworld = Traits :: dimensionworld };

    typedef typename Traits :: template Codim< codimension > :: Geometry Geometry;

    static const bool fake = false;

  private:
    typedef typename Traits :: HostGrid HostGrid;
    typedef typename Traits :: CoordFunction CoordFunction;

    template< class, bool > friend class GeometryGridEntityPointer;
    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class, class > friend class GeometryGridIdSet;
    template< class, class > friend class GeometryGridCommDataHandle;

  public:
    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

  private:
    typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

    typedef GeometryGridCoordVector< mydimension, Grid, fake > CoordVector;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

    const Grid *grid_;
    const HostEntity *hostEntity_;
    mutable MakeableGeometry geo_;

  public:
    GeometryGridEntityImpl ()
      : geo_( GeometryImpl() )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( GeometryImpl() )
    {}

  public:
    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl & );

  public:
    GeometryType type () const
    {
      return hostEntity().type();
    }

    int level () const
    {
      return hostEntity().level();
    }

    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

    const Geometry &geometry () const
    {
      GeometryImpl &geo = Grid :: getRealImplementation( geo_ );
      if( !geo )
      {
        CoordVector coords( hostEntity.geometry(), grid().coordFunction() );
        geo = GeometryImpl( type(), coords );
      }
      return geo_;
    }

    const HostEntity &hostEntity () const
    {
      return *hostEntity_;
    }

    const Grid &grid () const
    {
      return *grid_;
    }

    void initialize ( const Grid &grid, const HostEntity &hostEntity )
    {
      grid_ = &grid;
      hostEntity_ = &hostEntity;
      Grid :: getRealImplementation( geo_ ) = GeometryImpl();
    }

  private:
    template< class IndexSet >
    typename IndexSet :: IndexType index ( const IndexSet &indexSet ) const
    {
      return indexSet.template index< codimension >( hostEntity() );
    }

    template< class IdSet >
    typename IdSet :: IdType id ( const IdSet &idSet ) const
    {
      return idSet.template id< codimension >( hostEntity() );
    }
  };



  // GeometryGridEntityImpl (fake)
  // -----------------------------

  template< int codim, class Grid >
  class GeometryGridEntityImpl< codim, Grid, true >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

  public:
    typedef typename Traits :: ctype ctype;

    enum { codimension = codim };
    enum { dimension = Traits :: dimension };
    enum { mydimension = dimension - codimension };
    enum { dimensionworld = Traits :: dimensionworld };

    typedef typename Traits :: template Codim< codimension > :: Geometry Geometry;

    static const bool fake = true;

  private:
    typedef typename Traits :: HostGrid HostGrid;
    typedef typename Traits :: CoordFunction CoordFunction;

    template< class, bool > friend class GeometryGridEntityPointer;
    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class, class > friend class GeometryGridIdSet;
    template< class, class > friend class GeometryGridCommDataHandle;

  public:
    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;

    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

  private:
    typedef typename HostGrid :: template Codim< 0 > :: Geometry HostGeometry;
    typedef typename HostGrid :: template Codim< dimension > :: EntityPointer
    HostVertexPointer;

    typedef GeometryGridCoordVector< mydimension, Grid, fake > CoordVector;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

    const Grid *grid_;
    const HostElement *hostElement_;
    unsigned int subEntity_;
    mutable Geometry geo_;

  public:
    GeometryGridEntityImpl ()
      : geo_( GeometryImpl() )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostElement_( other.hostElement_ ),
        subEntity_( other.subEntity_ ),
        geo_( GeometryImpl() )
    {}

  private:
    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl & );

  public:
    GeometryType type () const
    {
      const ReferenceElement< ctype, dimension > &refElement
        = ReferenceElements< ctype, dimension > :: general( hostElement().type() );
      return refElement.type( subEntity_, codimension );
    }

    int level () const
    {
      return hostElement().level();
    }

    PartitionType partitionType () const
    {
      if( !(Capabilities :: isParallel< HostGrid > :: v) )
        return InteriorEntity;

      const ReferenceElement< ctype, dimension > &refElement
        = ReferenceElements< ctype, dimension > :: general( hostElement().type() );

      PartitionType type = vertexPartitionType( refElement, 0 );
      if( (type == InteriorEntity) || (type == OverlapEntity)
          || (type == GhostEntity) )
        return type;

      const int numVertices = refElement.size( subEntity_, codimension, dimension );
      for( int i = 1; i < numVertices; ++i )
      {
        PartitionType vtxType = vertexPartitionType( refElement, i );
        if( (vtxType == InteriorEntity) || (vtxType == OverlapEntity)
            || (vtxType == GhostEntity) )
          return vtxType;
        assert( type == vtxType );
      }
      assert( (type == BorderEntity) || (type == FrontEntity) );
      return type;
    }

    const Geometry &geometry () const
    {
      GeometryImpl &geo = Grid :: getRealImplementation( geo_ );
      if( !geo )
      {
        CoordVector coords( hostEntity.geometry(), subEntity_, grid().coordFunction() );
        geo = GeometryImpl( type(), coords );
      }
      return geo_;
    }

    const HostEntity &hostEntity () const
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

    const Grid &grid () const
    {
      return *grid_;
    }

    void initialize ( const Grid &grid, const HostEntity &hostEntity )
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

    void initialize ( const Grid &grid, const HostElement &hostElement, int subEntity )
    {
      grid_ = &grid;
      hostElement_ = &hostElement;
      subEntity_ = subEntity;
      Grid :: getRealImplementation( geo_ ) = GeometryImpl();
    }

  private:
    PartitionType
    vertexPartitionType ( const ReferenceElement< ctype, dimension > &refElement,
                          int i ) const
    {
      const int j = refElement.subEntity( subEntity_, codimension, 0, dimension );
      return hostElement().template entity< dimension >( j )->partitionType();
    }

    const HostElement &hostElement () const
    {
      return *hostElement_;
    }

    template< class IndexSet >
    typename IndexSet :: IndexType index ( const IndexSet &indexSet ) const
    {
      return indexSet.template subIndex< codimension >( hostElement(), subEntity_ );
    }

    template< class IdSet >
    typename IdSet :: IdType id ( const IdSet &idSet ) const
    {
      return idSet.template subId< codimension >( hostElement(), subEntity_ );
    }
  };



  // GeometryGridEntityImpl for codimension 0
  // ----------------------------------------

  template< class Grid >
  class GeometryGridEntityImpl< 0, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

  public:
    typedef typename Traits :: ctype ctype;

    enum { codimension = 0 };
    enum { dimension = Traits :: dimension };
    enum { mydimension = dimension - codimension };
    enum { dimensionworld = Traits :: dimensionworld };

    typedef typename Traits :: template Codim< codimension > :: EntityPointer EntityPointer;
    typedef typename Traits :: template Codim< codimension > :: Geometry Geometry;
    typedef typename Traits :: template Codim< codimension > :: LocalGeometry LocalGeometry;

    typedef typename Traits :: HierarchicIterator HierarchicIterator;
    typedef typename Traits :: LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Traits :: LevelIntersectionIterator LevelIntersectionIterator;

    static const bool fake = false;

  private:
    typedef typename Traits :: HostGrid HostGrid;
    typedef typename Traits :: CoordFunction CoordFunction;

    template< class, bool > friend class GeometryGridEntityPointer;
    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class, class > friend class GeometryGridIdSet;
    template< class, class > friend class GeometryGridCommDataHandle;

  public:
    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

  private:
    typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

    typedef GeometryGridCoordVector< mydimension, Grid, fake > CoordVector;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    const HostEntity *hostEntity_;
    mutable MakeableGeometry geo_;

  public:
    GeometryGridEntityImpl ()
      : geo_( GeometryImpl() )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( GeometryImpl() )
    {}

  private:
    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl & );

  public:
    GeometryType type () const
    {
      return hostEntity().type();
    }

    int level () const
    {
      return hostEntity().level();
    }

    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

    const Geometry &geometry () const
    {
      GeometryImpl &geo = Grid :: getRealImplementation( geo_ );

      if( !geo )
      {
        CoordVector coords( hostEntity().geometry(), grid().coordFunction() );
        geo = GeometryImpl( type(), coords );
      }
      return geo_;
    }

    template< int codim >
    int count () const
    {
      return hostEntity().template count< codim >();
    }

    template< int codim >
    typename Grid :: template Codim< codim > :: EntityPointer
    entity ( int i ) const
    {
      typedef typename Grid :: template Codim< codim > :: EntityPointer EntityPointer;
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;

      EntityPointerImpl impl( *grid_, hostEntity(), i );
      return MakeableEntityPointer( impl );
    }

    LevelIntersectionIterator ilevelbegin () const
    {
      typedef MakeableInterfaceObject< LevelIntersectionIterator >
      MakeableLevelIntersectionIterator;
      typedef typename MakeableLevelIntersectionIterator :: ImplementationType
      LevelIntersectionIteratorImpl;
      LevelIntersectionIteratorImpl impl( *this, hostEntity().ilevelbegin() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LevelIntersectionIterator ilevelend () const
    {
      typedef MakeableInterfaceObject< LevelIntersectionIterator >
      MakeableLevelIntersectionIterator;
      typedef typename MakeableLevelIntersectionIterator :: ImplementationType
      LevelIntersectionIteratorImpl;
      LevelIntersectionIteratorImpl impl( *this, hostEntity().ilevelend() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafbegin () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( *this, hostEntity().ileafbegin() );
      return MakeableLeafIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafend () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( *this, hostEntity().ileafend() );
      return MakeableLeafIntersectionIterator( impl );
    }

    bool hasBoundaryIntersections () const
    {
      return hostEntity().hasBoundaryIntersections();
    }

    bool isLeaf () const
    {
      return hostEntity().isLeaf();
    }

    EntityPointer father () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostEntity().father() ) );
    }

    const LocalGeometry &geometryInFather () const
    {
      return hostEntity().geometryInFather();
    }

    HierarchicIterator hbegin ( int maxLevel ) const
    {
      typedef MakeableInterfaceObject< HierarchicIterator > MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *grid_, hostEntity().hbegin( maxLevel ) );
      return MakeableIterator( impl );
    }

    HierarchicIterator hend ( int maxLevel ) const
    {
      typedef MakeableInterfaceObject< HierarchicIterator > MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *grid_, hostEntity().hend( maxLevel ) );
      return MakeableIterator( impl );
    }

    bool isRegular () const
    {
      return hostEntity().isRegular();
    }

    bool wasRefined () const
    {
      return hostEntity().wasRefined();
    }

    bool mightBeCoarsened () const
    {
      return hostEntity().mightBeCoarsened();
    }

    const HostEntity &hostEntity () const
    {
      return *hostEntity_;
    }

    const Grid &grid () const
    {
      return *grid_;
    }

    void initialize ( const Grid &grid, const HostEntity &hostEntity )
    {
      grid_ = &grid;
      hostEntity_ = &hostEntity;
      Grid :: getRealImplementation( geo_ ) = GeometryImpl();
    }

  private:
    template< class IndexSet >
    typename IndexSet :: IndexType index ( const IndexSet &indexSet ) const
    {
      return indexSet.template index< codimension >( hostEntity() );
    }

    template< class IdSet >
    typename IdSet :: IdType id ( const IdSet &idSet ) const
    {
      return idSet.template id< codimension >( hostEntity() );
    }
  };



  template< int codim, int dim, class Grid >
  class GeometryGridEntity
    : public GeometryGridEntityImpl< codim, Grid >
  {
    typedef GeometryGridEntityImpl< codim, Grid > Base;

  public:
    GeometryGridEntity ()
      : Base()
    {}

    GeometryGridEntity ( const Grid &grid )
      : Base( grid )
    {}
  };



  template< class Entity >
  class GeometryGridEntityWrapper;

  template< int codim, int dim, class Grid >
  class GeometryGridEntityWrapper< Entity< codim, dim, Grid, GeometryGridEntity > >
    : public Entity< codim, dim, Grid, GeometryGridEntity >
  {
    typedef Entity< codim, dim, Grid, GeometryGridEntity > Base;

  protected:
    using Base :: getRealImp;

  public:
    typedef GeometryGridEntity< codim, dim, Grid > Implementation;

    typedef typename Implementation :: HostEntity HostEntity;
    typedef typename Implementation :: HostElement HostElement;

    GeometryGridEntityWrapper ()
      : Base( Implementation() )
    {}

    void initialize ( const Grid &grid, const HostEntity &hostEntity )
    {
      getRealImp().initialize( grid, hostEntity );
    }

    void initialize ( const Grid &grid, const HostElement &hostElement, int subEntity )
    {
      getRealImp().initialize( grid, hostElement, subEntity );
    }
  };

}

#endif
