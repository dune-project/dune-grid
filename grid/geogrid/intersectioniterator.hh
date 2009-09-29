// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <vector>

#include <dune/grid/geogrid/entitypointer.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< class Grid, class HostIntersection >
  class GeometryGridIntersection;

  template< class Grid >
  class GeometryGridLeafIntersection;

  template< class Grid >
  class GeometryGridLevelIntersection;

  template< class Grid >
  class GeometryGridLeafIntersectionIterator;

  template< class Grid >
  class GeometryGridLevelIntersectionIterator;



  // GeometryGridIntersectionIterator
  // --------------------------------

  template< class HostGrid, class CoordFunction, class HostIntersection >
  class GeometryGridIntersection
  < const GeometryGrid< HostGrid, CoordFunction >, HostIntersection >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostIntersection :: Geometry HostGeometry;
    typedef typename HostIntersection :: LocalGeometry HostLocalGeometry;

  public:
    typedef typename Grid :: ctype ctype;

    enum { dimension = Grid :: dimension };
    enum { dimensionworld = Grid :: dimensionworld };

    typedef typename Grid :: template Codim< 0 > :: Entity Entity;
    typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;
    typedef typename Grid :: template Codim< 1 > :: Geometry Geometry;
    typedef typename Grid :: template Codim< 1 > :: LocalGeometry LocalGeometry;

  private:
    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    const HostIntersection *hostIntersection_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;
    mutable LocalGeometry *geoSelf_;
    mutable LocalGeometry *geoNeighbor_;

  public:
    GeometryGridIntersection ( const Grid &grid )
      : grid_( &grid ),
        hostIntersection_( 0 ),
        geo_( 0 ),
        geoSelf_( 0 ),
        geoNeighbor_( 0 )
    {}

    GeometryGridIntersection ( const GeometryGridIntersection &other )
      : grid_( other.grid_ ),
        hostIntersection_( other.hostIntersection_ ),
        geo_( 0 ),
        geoSelf_( 0 ),
        geoNeighbor_( 0 )
    {}

    ~GeometryGridIntersection ()
    {
      if( geo_ != 0 )
        delete geo_;
      if( geoSelf_ != 0 )
        delete geoSelf_;
      if( geoNeighbor_ != 0 )
        delete geoNeighbor_;
    }

    EntityPointer inside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIntersection().inside() ) );
    }

    EntityPointer outside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIntersection().outside() ) );
    }

    bool boundary () const
    {
      return hostIntersection().boundary ();
    }

    bool neighbor () const
    {
      return hostIntersection().neighbor();
    }

    int boundaryId () const
    {
      return hostIntersection().boundaryId();
    }

    const LocalGeometry &intersectionSelfLocal () const
    {
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry::ImplementationType LocalGeometryImpl;

      if( geoSelf_ == 0 )
      {
        LocalGeometryImpl impl( hostIntersection().intersectionSelfLocal() );
        geoSelf_ = new MakeableLocalGeometry( impl );
      }
      return *geoSelf_;
    }

    const LocalGeometry &intersectionNeighborLocal () const
    {
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry::ImplementationType LocalGeometryImpl;

      if( geoNeighbor_ == 0 )
      {
        LocalGeometryImpl impl( hostIntersection().intersectionNeighborLocal() );
        geoNeighbor_ = new MakeableLocalGeometry( impl );
      }
      return *geoNeighbor_;
    }

    const Geometry &intersectionGlobal () const
    {
      if( geo_ == 0 )
      {
        const HostGeometry &hostGeo = hostIntersection().intersectionGlobal();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( hostGeo.type(), corners_ ) );
      }
      return *geo_;
    }

    int numberInSelf () const
    {
      return hostIntersection().numberInSelf();
    }

    int numberInNeighbor () const
    {
      return hostIntersection().numberInNeighbor();
    }

    FieldVector< ctype, dimensionworld >
    integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      typedef typename Grid :: template Codim< 0 > :: Geometry Geometry;
      EntityPointer insideEntity = inside();
      const Geometry &geo = insideEntity->geometry();
      FieldVector< ctype, dimension > x( intersectionSelfLocal().global( local ) );
      return Grid :: getRealImplementation( geo ).normal( numberInSelf(), x );
    }

    FieldVector< ctype, dimensionworld >
    outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return integrationOuterNormal( local );
    }

    FieldVector< ctype, dimensionworld >
    unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      FieldVector< ctype, dimensionworld > normal = outerNormal( local );
      normal *= (ctype( 1 ) / normal.two_norm());
      return normal;
    }

  protected:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }

    bool isValid () const
    {
      return (hostIntersection_ != 0);
    }

    const HostIntersection &hostIntersection () const
    {
      assert( isValid() );
      return *hostIntersection_;
    }

    void invalidate ()
    {
      hostIntersection_ = 0;
    }

    void setToTarget ( const HostIntersection &hostIntersection )
    {
      hostIntersection_ = &hostIntersection;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      if( geoSelf_ != 0 )
      {
        delete geoSelf_;
        geoSelf_ = 0;
      }

      if( geoNeighbor_ != 0 )
      {
        delete geoNeighbor_;
        geoNeighbor_ = 0;
      }
    }
  };


  // GeometryGridLeafIntersection
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLeafIntersection< const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridIntersection
      < const GeometryGrid< HostGrid, CoordFunction >,
          typename HostGrid :: Traits :: LeafIntersection >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    friend class GeometryGridLeafIntersectionIterator< const Grid >;

    typedef typename HostGrid :: Traits :: LeafIntersection HostIntersection;

    typedef GeometryGridIntersection< const Grid, HostIntersection > Base;

  public:
    GeometryGridLeafIntersection ( const Grid &grid )
      : Base( grid )
    {}
  };



  // GeometryGridLevelIntersection
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLevelIntersection< const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridIntersection
      < const GeometryGrid< HostGrid, CoordFunction >,
          typename HostGrid :: Traits :: LevelIntersection >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    friend class GeometryGridLevelIntersectionIterator< const Grid >;

    typedef typename HostGrid :: Traits :: LevelIntersection HostIntersection;

    typedef GeometryGridIntersection< const Grid, HostIntersection > Base;

  public:
    GeometryGridLevelIntersection ( const Grid &grid )
      : Base( grid )
    {}
  };



  // GeometryGridLeafIntersectionIterator
  // ------------------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLeafIntersectionIterator
  < const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: Traits :: LeafIntersectionIterator
    HostIntersectionIterator;

  public:
    typedef typename Grid :: Traits :: LeafIntersection Intersection;

  private:
    typedef MakeableInterfaceObject< Intersection > MakeableIntersection;
    typedef typename MakeableIntersection :: ImplementationType IntersectionImpl;

    HostIntersectionIterator hostIterator_;
    mutable MakeableIntersection intersection_;

  public:
    GeometryGridLeafIntersectionIterator
      ( const Grid &grid,
      const HostIntersectionIterator &hostIterator )
      : hostIterator_( hostIterator ),
        intersection_( IntersectionImpl( grid ) )
    {}

    bool equals ( const GeometryGridLeafIntersectionIterator &other ) const
    {
      return (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;
      Grid :: getRealImplementation( intersection_ ).invalidate();
    }

    const Intersection &dereference () const
    {
      IntersectionImpl &impl = Grid :: getRealImplementation( intersection_ );
      if( !impl.isValid() )
        impl.setToTarget( *hostIterator_ );
      return intersection_;
    }
  };



  // GeometryGridLevelIntersectionIterator
  // -------------------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLevelIntersectionIterator
  < const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: Traits :: LevelIntersectionIterator
    HostIntersectionIterator;

  public:
    typedef typename Grid :: Traits :: LevelIntersection Intersection;

  private:
    typedef MakeableInterfaceObject< Intersection > MakeableIntersection;
    typedef typename MakeableIntersection :: ImplementationType IntersectionImpl;

    HostIntersectionIterator hostIterator_;
    mutable MakeableIntersection intersection_;

  public:
    GeometryGridLevelIntersectionIterator
      ( const Grid &grid,
      const HostIntersectionIterator &hostIterator )
      : hostIterator_( hostIterator ),
        intersection_( IntersectionImpl( grid ) )
    {}

    bool equals ( const GeometryGridLevelIntersectionIterator &other ) const
    {
      return (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;
      Grid :: getRealImplementation( intersection_ ).invalidate();
    }

    const Intersection &dereference () const
    {
      IntersectionImpl &impl = Grid :: getRealImplementation( intersection_ );
      if( !impl.isValid() )
        impl.setToTarget( *hostIterator_ );
      return intersection_;
    }
  };

}

#endif
