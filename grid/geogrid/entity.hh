// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/referenceelements.hh>

#include <dune/grid/geogrid/capabilities.hh>

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

    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
    typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    const HostEntity *hostEntity_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;

  public:
    explicit GeometryGridEntityImpl ( const Grid &grid )
      : grid_( &grid ),
        hostEntity_( 0 ),
        geo_( 0 )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( 0 )
    {}

    ~GeometryGridEntityImpl ()
    {
      if( geo_ != 0 )
        delete geo_;
    }

    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl &other )
    {
      if( this == &other )
        return *this;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      grid_ = other.grid_;
      hostEntity_ = other.hostEntity_;
      return *this;
    }

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
      if( geo_ == 0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        const CoordFunction &coordFunction = grid().coordFunction();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction.evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
    }

    const HostEntity &hostEntity () const
    {
      assert( isValid() );
      return *hostEntity_;
    }

    const Grid &grid () const
    {
      return *grid_;
    }

  private:
    bool isValid () const
    {
      return (hostEntity_ != 0);
    }

    void invalidate ()
    {
      hostEntity_ = 0;
    }

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

    void setToTarget( const HostEntity &hostEntity )
    {
      hostEntity_ = &hostEntity;
      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }
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

    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;

    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;
    typedef typename HostGrid :: template Codim< 0 > :: Geometry HostGeometry;
    typedef typename HostGrid :: template Codim< dimension > :: EntityPointer
    HostVertexPointer;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    const HostElement *hostElement_;
    unsigned int subEntity_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;

  public:
    explicit GeometryGridEntityImpl ( const Grid &grid )
      : grid_( &grid ),
        hostElement_( 0 ),
        geo_( 0 )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostElement_( other.hostElement_ ),
        subEntity_( other.subEntity_ ),
        geo_( 0 )
    {}

    ~GeometryGridEntityImpl ()
    {
      if( geo_ != 0 )
        delete geo_;
    }

    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl &other )
    {
      if( this == &other )
        return *this;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      grid_ = other.grid_;
      hostElement_ = other.hostElement_;
      subEntity_ = other.subEntity_;
      return *this;
    }

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
      if( geo_ == 0 )
      {
        const ReferenceElement< ctype, dimension > &refElement
          = ReferenceElements< ctype, dimension > :: general( hostElement().type() );
        corners_.resize( refElement.size( subEntity_, codimension, dimension ) );

        const HostGeometry &hostGeo = hostElement().geometry();
        const CoordFunction &coordFunction = grid().coordFunction();
        for( unsigned int i = 0; i < corners_.size(); ++i )
        {
          const int j = refElement.subEntity( subEntity_, codimension, i, dimension );
          coordFunction.evaluate( hostGeo[ j ], corners_[ i ] );
        }
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
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
      assert( isValid() );
      return *hostElement_;
    }

    bool isValid () const
    {
      return (hostElement_ != 0);
    }

    void invalidate ()
    {
      hostElement_ = 0;
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

    void setToTarget ( const HostEntity &hostEntity )
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

    void setToTarget ( const HostElement &hostElement, int subEntity )
    {
      hostElement_ = &hostElement;
      subEntity_ = subEntity;
      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }
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

    typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
    typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    const HostEntity *hostEntity_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;
    mutable LocalGeometry *geoInFather_;

  public:
    explicit GeometryGridEntityImpl ( const Grid &grid )
      : grid_( &grid ),
        hostEntity_( 0 ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    GeometryGridEntityImpl ( const GeometryGridEntityImpl &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    ~GeometryGridEntityImpl ()
    {
      if( geo_ != 0 )
        delete geo_;
      if( geoInFather_ != 0 )
        delete geoInFather_;
    }


    GeometryGridEntityImpl &operator= ( const GeometryGridEntityImpl &other )
    {
      if( this == &other )
        return *this;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      if( geoInFather_ != 0 )
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }

      grid_ = other.grid_;
      hostEntity_ = other.hostEntity_;
      return *this;
    }

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
      if( geo_ == 0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        const CoordFunction &coordFunction = grid().coordFunction();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction.evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
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
      LevelIntersectionIteratorImpl impl( *grid_, hostEntity().ilevelbegin() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LevelIntersectionIterator ilevelend () const
    {
      typedef MakeableInterfaceObject< LevelIntersectionIterator >
      MakeableLevelIntersectionIterator;
      typedef typename MakeableLevelIntersectionIterator :: ImplementationType
      LevelIntersectionIteratorImpl;
      LevelIntersectionIteratorImpl impl( *grid_, hostEntity().ilevelend() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafbegin () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( *grid_, hostEntity().ileafbegin() );
      return MakeableLeafIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafend () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( *grid_, hostEntity().ileafend() );
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
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry :: ImplementationType LocalGeometryImpl;

      if( geoInFather_ == 0 )
        geoInFather_ = new MakeableLocalGeometry( LocalGeometryImpl( hostEntity().geometryInFather() ) );
      return *geoInFather_;
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

  private:
    bool isValid () const
    {
      return (hostEntity_ != 0);
    }

    void invalidate ()
    {
      hostEntity_ = 0;
    }

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

    void setToTarget( const HostEntity &hostEntity )
    {
      hostEntity_ = &hostEntity;
      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }
      if( geoInFather_ != 0 )
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
    }
  };



  template< int codim, int dim, class Grid >
  class GeometryGridEntity
    : public GeometryGridEntityImpl< codim, Grid >
  {
    typedef GeometryGridEntityImpl< codim, Grid > Base;

  public:
    GeometryGridEntity ( const Grid &grid )
      : Base( grid )
    {}
  };



  template< int codim, int dim, class Grid >
  class GeometryGridEntityWrapper
    : public Entity< codim, dim, Grid, GeometryGridEntity >
  {
    typedef Entity< codim, dim, Grid, GeometryGridEntity > Base;

  public:
    typedef GeometryGridEntity< codim, dim, Grid > Implementation;

    GeometryGridEntityWrapper ( const Grid &grid )
      : Base( Implementation( grid ) )
    {}
  };

}

#endif
