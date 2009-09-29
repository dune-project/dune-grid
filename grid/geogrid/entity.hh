// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/referenceelements.hh>

#include <dune/grid/geogrid/capabilities.hh>
#include <dune/grid/geogrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< int codim, class Grid,
        bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v) >
    class EntityImpl;

    template< int codim, int dim, class Grid >
    class Entity;



    // EntityImpl (real)
    // -----------------

    template< int codim, class Grid >
    class EntityImpl< codim, Grid, false >
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

    public:
      typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
      typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
      typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

    private:
      typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

      typedef GeoGrid :: CoordVector< mydimension, Grid, fake > CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

      const Grid *grid_;
      const HostEntity *hostEntity_;
      mutable MakeableGeometry geo_;

    public:
      EntityImpl ()
        : geo_( GeometryImpl() )
      {}

      EntityImpl ( const EntityImpl &other )
        : grid_( other.grid_ ),
          hostEntity_( other.hostEntity_ ),
          geo_( GeometryImpl() )
      {}

    public:
      EntityImpl &operator= ( const EntityImpl & );

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


      // Implementation Stuff

      const Grid &grid () const
      {
        return *grid_;
      }

      const HostEntity &hostEntity () const
      {
        return *hostEntity_;
      }

      void initialize ( const Grid &grid, const HostEntity &hostEntity )
      {
        grid_ = &grid;
        hostEntity_ = &hostEntity;
        Grid :: getRealImplementation( geo_ ) = GeometryImpl();
      }

      template< class HostIndexSet >
      typename HostIndexSet :: IndexType
      index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.template index< codimension >( hostEntity() );
      }

      template< class HostIdSet >
      typename HostIdSet :: IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.template id< codimension >( hostEntity() );
      }
    };



    // EntityImpl (fake)
    // -----------------

    template< int codim, class Grid >
    class EntityImpl< codim, Grid, true >
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

    public:
      typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
      typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;

      typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

    private:
      typedef typename HostGrid :: template Codim< 0 > :: Geometry HostGeometry;
      typedef typename HostGrid :: template Codim< dimension > :: EntityPointer
      HostVertexPointer;

      typedef GeoGrid :: CoordVector< mydimension, Grid, fake > CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

      const Grid *grid_;
      const HostElement *hostElement_;
      unsigned int subEntity_;
      mutable Geometry geo_;

    public:
      EntityImpl ()
        : geo_( GeometryImpl() )
      {}

      EntityImpl ( const EntityImpl &other )
        : grid_( other.grid_ ),
          hostElement_( other.hostElement_ ),
          subEntity_( other.subEntity_ ),
          geo_( GeometryImpl() )
      {}

    private:
      EntityImpl &operator= ( const EntityImpl & );

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
          CoordVector coords( hostElement().geometry(), subEntity_, grid().coordFunction() );
          geo = GeometryImpl( type(), coords );
        }
        return geo_;
      }


      // Implementation Stuff

      const Grid &grid () const
      {
        return *grid_;
      }

      const HostElement &hostElement () const
      {
        return *hostElement_;
      }

      const HostEntity &hostEntity () const
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

      template< class HostIndexSet >
      typename HostIndexSet :: IndexType index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.template subIndex< codimension >( hostElement(), subEntity_ );
      }

      template< class HostIdSet >
      typename HostIdSet :: IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.template subId< codimension >( hostElement(), subEntity_ );
      }

    private:
      PartitionType
      vertexPartitionType ( const ReferenceElement< ctype, dimension > &refElement,
                            int i ) const
      {
        const int j = refElement.subEntity( subEntity_, codimension, 0, dimension );
        return hostElement().template entity< dimension >( j )->partitionType();
      }
    };



    // EntityImpl for codimension 0
    // ----------------------------

    template< class Grid >
    class EntityImpl< 0, Grid, false >
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

    public:
      typedef typename HostGrid :: template Codim< codimension > :: Entity HostEntity;
      typedef typename HostGrid :: template Codim< codimension > :: EntityPointer HostEntityPointer;
      typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

    private:
      typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

      typedef GeoGrid :: CoordVector< mydimension, Grid, fake > CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
      typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

      const Grid *grid_;
      const HostEntity *hostEntity_;
      mutable MakeableGeometry geo_;

    public:
      EntityImpl ()
        : geo_( GeometryImpl() )
      {}

      EntityImpl ( const EntityImpl &other )
        : grid_( other.grid_ ),
          hostEntity_( other.hostEntity_ ),
          geo_( GeometryImpl() )
      {}

    private:
      EntityImpl &operator= ( const EntityImpl & );

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


      // Implementation Stuff

      const Grid &grid () const
      {
        return *grid_;
      }

      const HostEntity &hostEntity () const
      {
        return *hostEntity_;
      }

      void initialize ( const Grid &grid, const HostEntity &hostEntity )
      {
        grid_ = &grid;
        hostEntity_ = &hostEntity;
        Grid :: getRealImplementation( geo_ ) = GeometryImpl();
      }

      template< class HostIndexSet >
      typename HostIndexSet :: IndexType index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.template index< codimension >( hostEntity() );
      }

      template< class HostIdSet >
      typename HostIdSet :: IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.template id< codimension >( hostEntity() );
      }
    };



    template< int codim, int dim, class Grid >
    class Entity
      : public EntityImpl< codim, Grid >
    {
      typedef EntityImpl< codim, Grid > Base;

    public:
      Entity ()
        : Base()
      {}

      Entity ( const Grid &grid )
        : Base( grid )
      {}
    };



    template< class Entity >
    class EntityWrapper;

    template< int codim, int dim, class Grid >
    class EntityWrapper< Dune :: Entity< codim, dim, Grid, Entity > >
      : public Dune :: Entity< codim, dim, Grid, Entity >
    {
      typedef Dune :: Entity< codim, dim, Grid, Entity > Base;

    protected:
      using Base :: getRealImp;

    public:
      typedef Entity< codim, dim, Grid > Implementation;

      typedef typename Implementation :: HostEntity HostEntity;
      typedef typename Implementation :: HostElement HostElement;

      EntityWrapper ()
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

}

#endif
