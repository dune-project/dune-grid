// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geogrid/entitypointer.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< class Grid >
  class GeometryGridLeafIntersectionIterator;

  template< class Grid >
  class GeometryGridLevelIntersectionIterator;



  // GeometryGridLeafIntersectionIterator
  // ------------------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLeafIntersectionIterator< const GeometryGrid< HostGrid, CoordFunction > >
    : public IntersectionIteratorDefaultImplementation
      < const GeometryGrid< HostGrid, CoordFunction >, GeometryGridLeafIntersectionIterator >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename Grid :: ctype ctype;

    enum { dimension = Grid :: dimension };
    enum { dimensionworld = Grid :: dimensionworld };

    typedef typename HostGrid :: Traits :: LeafIntersectionIterator
    HostLeafIntersectionIterator;
    typedef typename HostGrid :: template Codim< 1 > :: Geometry HostGeometry;

  public:
    typedef typename Grid :: template Codim< 0 > :: Entity Entity;
    typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;
    typedef typename Grid :: template Codim< 1 > :: Geometry Geometry;
    typedef typename Grid :: template Codim< 1 > :: LocalGeometry LocalGeometry;

  private:
    typedef MakeableInterfaceObject<Geometry> MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    HostLeafIntersectionIterator hostIterator_;
    mutable std::vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;
    mutable LocalGeometry *geoSelf_;
    mutable LocalGeometry *geoNeighbor_;

  public:
    GeometryGridLeafIntersectionIterator
      ( const Grid *grid,
      const HostLeafIntersectionIterator &hostIterator )
      : grid_( grid ),
        hostIterator_( hostIterator ),
        geo_( 0 ),
        geoSelf_( 0 ),
        geoNeighbor_( 0 )
    {}

    ~GeometryGridLeafIntersectionIterator ()
    {}

    bool equals ( const GeometryGridLeafIntersectionIterator &other ) const
    {
      return (grid_ == other.grid_) && (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;

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

    EntityPointer inside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIterator_->inside() ) );
    }

    EntityPointer outside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIterator_->outside() ) );
    }

    bool boundary () const
    {
      return hostIterator_->boundary ();
    }

    bool neighbor () const
    {
      return hostIterator_->neighbor();
    }

    int boundaryId () const
    {
      return hostIterator_->boundaryId();
    }

    const LocalGeometry &intersectionSelfLocal () const
    {
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry::ImplementationType LocalGeometryImpl;

      if( geoSelf_ == 0 )
      {
        LocalGeometryImpl impl( hostIterator_->intersectionSelfLocal() );
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
        LocalGeometryImpl impl( hostIterator_->intersectionNeighborLocal() );
        geoNeighbor_ = new MakeableLocalGeometry( impl );
      }
      return *geoNeighbor_;
    }

    const Geometry &intersectionGlobal () const
    {
      if( geo_ == 0 )
      {
        const HostGeometry &hostGeo = hostIterator_->intersectionGlobal();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( hostGeo.type(), corners_ ) );
      }
      return *geo_;
    }

    int numberInSelf () const
    {
      return hostIterator_->numberInSelf();
    }

    int numberInNeighbor () const
    {
      return hostIterator_->numberInNeighbor();
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

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }
  };



  // GeometryGridLevelIntersectionIterator
  // -------------------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLevelIntersectionIterator< const GeometryGrid< HostGrid, CoordFunction > >
    : public IntersectionIteratorDefaultImplementation
      < const GeometryGrid< HostGrid, CoordFunction >, GeometryGridLevelIntersectionIterator >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename Grid :: ctype ctype;

    enum { dimension = Grid :: dimension };
    enum { dimensionworld = Grid :: dimensionworld };

    typedef typename HostGrid :: Traits :: LevelIntersectionIterator
    HostLevelIntersectionIterator;
    typedef typename HostGrid :: template Codim< 1 > :: Geometry HostGeometry;

  public:
    typedef typename Grid :: template Codim< 0 > :: Entity Entity;
    typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;
    typedef typename Grid :: template Codim< 1 > :: Geometry Geometry;
    typedef typename Grid :: template Codim< 1 > :: LocalGeometry LocalGeometry;

  private:
    typedef MakeableInterfaceObject<Geometry> MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const Grid *grid_;
    HostLevelIntersectionIterator hostIterator_;
    mutable std::vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;
    mutable LocalGeometry *geoSelf_;
    mutable LocalGeometry *geoNeighbor_;

  public:
    GeometryGridLevelIntersectionIterator
      ( const Grid *grid,
      const HostLevelIntersectionIterator &hostIterator )
      : grid_( grid ),
        hostIterator_( hostIterator ),
        geo_( 0 ),
        geoSelf_( 0 ),
        geoNeighbor_( 0 )
    {}

    ~GeometryGridLevelIntersectionIterator ()
    {}

    bool equals ( const GeometryGridLevelIntersectionIterator &other ) const
    {
      return (grid_ == other.grid_) && (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;

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

    EntityPointer inside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIterator_->inside() ) );
    }

    EntityPointer outside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( *grid_, hostIterator_->outside() ) );
    }

    bool boundary () const
    {
      return hostIterator_->boundary ();
    }

    bool neighbor () const
    {
      return hostIterator_->neighbor();
    }

    int boundaryId () const
    {
      return hostIterator_->boundaryId();
    }

    const LocalGeometry &intersectionSelfLocal () const
    {
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry::ImplementationType LocalGeometryImpl;

      if( geoSelf_ == 0 )
      {
        LocalGeometryImpl impl( hostIterator_->intersectionSelfLocal() );
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
        LocalGeometryImpl impl( hostIterator_->intersectionNeighborLocal() );
        geoNeighbor_ = new MakeableLocalGeometry( impl );
      }
      return *geoNeighbor_;
    }

    const Geometry &intersectionGlobal () const
    {
      if( geo_ == 0 )
      {
        const HostGeometry &hostGeo = hostIterator_->intersectionGlobal();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( hostGeo.type(), corners_ ) );
      }
      return *geo_;
    }

    int numberInSelf () const
    {
      return hostIterator_->numberInSelf();
    }

    int numberInNeighbor () const
    {
      return hostIterator_->numberInNeighbor();
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

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }
  };

}  // namespace Dune

#endif
