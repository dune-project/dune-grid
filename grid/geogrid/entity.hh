// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/referenceelements.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridFamily;

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;

  template<int codim, class GridImp>
  class GeometryGridEntityPointer;

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class GeometryGridLevelIterator;

  template<class GridImp>
  class GeometryGridLevelIntersectionIterator;

  template<class GridImp>
  class GeometryGridLeafIntersectionIterator;

  template<class GridImp>
  class GeometryGridHierarchicIterator;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class GeometryGridEntity;



  // GeometryGridEntity
  // ------------------

  template< int codim, int dim, class HostGrid, class CoordFunction >
  class GeometryGridEntity< codim, dim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef typename GeometryGridFamily< HostGrid, CoordFunction > :: Traits Traits;

    typedef typename Traits :: Grid Grid;

  public:
    typedef typename Traits :: ctype ctype;

    enum { codimension = codim };
    enum { dimension = dim };
    enum { mydimension = dimension - codimension };
    enum { dimensionworld = Traits :: dimensionworld };

    typedef typename Traits :: template Codim< codimension > :: Geometry Geometry;

  private:
    friend class GeometryGridEntityPointer< codimension, const Grid >;

    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class > friend class GeometryGridLocalIdSet;
    template< class > friend class GeometryGridGlobalIdSet;
    template< class, int > friend class IndexSetter;

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
    GeometryGridEntity( const Grid &grid )
      : grid_( &grid ),
        hostEntity_( 0 ),
        geo_( 0 )
    {}

    GeometryGridEntity( const GeometryGridEntity &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( 0 )
    {}

    ~GeometryGridEntity ()
    {
      if( geo_ != 0 )
        delete geo_;
    }

    GeometryGridEntity &operator= ( const GeometryGridEntity &other )
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
      if( geo_ ==0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
    }

    const HostEntity &hostEntity () const
    {
      assert( isValid() );
      return *hostEntity_;
    }

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }

    bool isValid () const
    {
      return (hostEntity_ != 0);
    }

    void invalidate ()
    {
      hostEntity_ = 0;
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


  // GeometryGridEntity for codimension 0
  // ------------------------------------

  template< int dim, class HostGrid, class CoordFunction >
  class GeometryGridEntity< 0, dim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef typename GeometryGridFamily< HostGrid, CoordFunction > :: Traits Traits;

    typedef typename Traits :: Grid Grid;

  public:
    typedef typename Traits :: ctype ctype;

    enum { codimension = 0 };
    enum { dimension = dim };
    enum { mydimension = dimension - codimension };
    enum { dimensionworld = Traits :: dimensionworld };

    typedef typename Traits :: template Codim< codimension > :: EntityPointer EntityPointer;
    typedef typename Traits :: template Codim< codimension > :: Geometry Geometry;
    typedef typename Traits :: template Codim< codimension > :: LocalGeometry LocalGeometry;

    typedef typename Traits :: HierarchicIterator HierarchicIterator;
    typedef typename Traits :: LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Traits :: LevelIntersectionIterator LevelIntersectionIterator;

  private:
    friend class GeometryGridEntityPointer< codimension, const Grid >;

    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class > friend class GeometryGridLocalIdSet;
    template< class > friend class GeometryGridGlobalIdSet;
    template< class, int > friend class IndexSetter;

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
    GeometryGridEntity ( const Grid &grid )
      : grid_( &grid ),
        hostEntity_( 0 ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    GeometryGridEntity( const GeometryGridEntity &other )
      : grid_( other.grid_ ),
        hostEntity_( other.hostEntity_ ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    ~GeometryGridEntity ()
    {
      if( geo_ != 0 )
        delete geo_;
      if( geoInFather_ != 0 )
        delete geoInFather_;
    }


    GeometryGridEntity &operator= ( const GeometryGridEntity &other )
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
      if( geo_ ==0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
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

      EntityPointerImpl impl( *grid_, hostEntity().template entity< codim >( i ) );
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
      return MakeableEntityPointer( EntityPointerImpl( grid_, hostEntity().father() ) );
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

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }

    bool isValid () const
    {
      return (hostEntity_ != 0);
    }

    void invalidate ()
    {
      hostEntity_ = 0;
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

}

#endif
