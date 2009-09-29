// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geogrid/entitypointer.hh>
//#include <dune/grid/geogrid/cache.hh>
#include <dune/grid/geogrid/cornerstorage.hh>
#include <dune/grid/geogrid/storage.hh>

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

  template< class Intersection >
  class GeometryGridIntersectionWrapper;

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
    typedef GeometryGridCoordVector< dimension-1, const Grid, false > CoordVector;

    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    // typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

    const EntityPointer *inside_;
    const HostIntersection *hostIntersection_;
    // mutable GeometryGridCache< GlobalCoordinate > corners_;
    mutable MakeableGeometry geo_;

  public:
    GeometryGridIntersection ()
      : geo_( GeometryImpl() )
    {}

    GeometryGridIntersection ( const GeometryGridIntersection &other )
      : inside_( other.inside_ ),
        hostIntersection_( other.hostIntersection_ ),
        geo_( GeometryImpl() )
    {}

    //EntityPointer inside () const
    const EntityPointer &inside () const
    {
      // typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      // typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      // return MakeableEntityPointer( EntityPointerImpl( grid(), hostIntersection().inside() ) );
      return *inside_;
    }

    EntityPointer outside () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( grid(), hostIntersection().outside() ) );
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
      return hostIntersection().intersectionSelfLocal();
    }

    const LocalGeometry &intersectionNeighborLocal () const
    {
      return hostIntersection().intersectionNeighborLocal();
    }

    const Geometry &intersectionGlobal () const
    {
      GeometryImpl &geo = Grid :: getRealImplementation( geo_ );
      if( !geo )
      {
        const HostGeometry &hostGeo = hostIntersection().intersectionGlobal();
        CoordVector coords( hostGeo, grid().coordFunction() );
        geo = GeometryImpl( hostGeo.type(), coords );

        /*
           const CoordFunction &coordFunction = grid().coordFunction();

           const unsigned int numCorners = hostGeo.corners();
           corners_.reserve( numCorners );
           for( unsigned int i = 0; i < numCorners; ++i )
           coordFunction.evaluate( hostGeo[ i ], corners_[ i ] );
           geo = GeometryImpl( hostGeo.type(), corners_ );
         */
      }
      return geo_;
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
      const Geometry &geo = inside()->geometry();
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

    void initialize( const EntityPointer &inside, const HostIntersection &hostIntersection )
    {
      inside_ = &inside;
      hostIntersection_ = &hostIntersection;
      Grid :: getRealImplementation( geo_ ) = GeometryImpl();
    }

  protected:
    const Grid &grid () const
    {
      return Grid :: getRealImplementation( inside() ).grid();
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
    typedef typename HostGrid :: Traits :: LeafIntersection HostIntersection;

    template< class > friend class GeometryGridIntersectionWrapper;
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
    typedef typename HostGrid :: Traits :: LevelIntersection HostIntersection;

    template< class > friend class GeometryGridIntersectionWrapper;
  };



  // GeomegryGridIntersectionWrapper
  // -------------------------------

  template< class Intersection >
  class GeometryGridIntersectionWrapper
    : public Intersection
  {
    typedef Intersection Base;

  protected:
    using Base :: getRealImp;

  public:
    typedef typename Intersection :: ImplementationType Implementation;

    typedef typename Implementation :: EntityPointer EntityPointer;
    typedef typename Implementation :: HostIntersection HostIntersection;

    GeometryGridIntersectionWrapper ()
      : Base( Implementation() )
    {}

    void initialize( const EntityPointer &inside,
                     const HostIntersection &hostIntersection )
    {
      getRealImp().initialize( inside, hostIntersection );
    }
  };



  // GeometryGridIntersectionIterator
  // --------------------------------

  template< class Traits >
  class GeometryGridIntersectionIterator
  {
    typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

  public:
    typedef typename Traits :: Intersection Intersection;
    typedef typename Traits :: GridTraits :: Grid Grid;

    typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;

  private:
    typedef GeometryGridIntersectionWrapper< Intersection > IntersectionWrapper;
    typedef GeometryGridStorage< IntersectionWrapper > IntersectionStorage;

    EntityPointer inside_;
    HostIntersectionIterator hostIterator_;
    mutable IntersectionWrapper *intersection_;

  public:
    template< class Entity >
    GeometryGridIntersectionIterator ( const Entity &inside,
                                       const HostIntersectionIterator &hostIterator )
      : inside_( inside.template entity< 0 >( 0 ) ),
        hostIterator_( hostIterator ),
        intersection_( 0 )
    {}

    GeometryGridIntersectionIterator  ( const GeometryGridIntersectionIterator &other )
      : inside_( other.inside_ ),
        hostIterator_( other.hostIterator_ ),
        intersection_( 0 )
    {}

    GeometryGridIntersectionIterator &
    operator= ( const GeometryGridIntersectionIterator &other )
    {
      inside_ = other.inside_;
      hostIterator_ = other.hostIterator_;
      update();
      return *this;
    }

    ~GeometryGridIntersectionIterator ()
    {
      IntersectionStorage :: free( intersection_ );
    }

    bool equals ( const GeometryGridIntersectionIterator &other ) const
    {
      return (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;
      update();
    }

    const Intersection &dereference () const
    {
      if( intersection_ == 0 )
      {
        intersection_ = IntersectionStorage :: alloc();
        intersection_->initialize( inside_, *hostIterator_ );
      }
      return *intersection_;
    }

  private:
    const Grid &grid () const
    {
      return Grid :: getRealImplementation( inside_ ).grid();
    }

    void update ()
    {
      IntersectionStorage :: free( intersection_ );
      intersection_ = 0;
    }
  };



  // GeometryGridLeafIntersectionIteratorTraits
  // ------------------------------------------

  template< class Grid >
  struct GeometryGridLeafIntersectionIteratorTraits
  {
    typedef typename remove_const< Grid > :: type :: Traits GridTraits;

    typedef typename GridTraits :: LeafIntersection Intersection;

    typedef typename GridTraits :: HostGrid :: Traits :: LeafIntersectionIterator
    HostIntersectionIterator;
  };



  // GeometryGridLeafIntersectionIterator
  // ------------------------------------

  template< class Grid >
  class GeometryGridLeafIntersectionIterator
    : public GeometryGridIntersectionIterator
      < GeometryGridLeafIntersectionIteratorTraits< Grid > >
  {
    typedef GeometryGridLeafIntersectionIteratorTraits< Grid > Traits;
    typedef GeometryGridIntersectionIterator< Traits > Base;

    typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

  public:
    typedef typename Traits :: Intersection Intersection;

  public:
    template< class Entity >
    GeometryGridLeafIntersectionIterator
      ( const Entity &inside,
      const HostIntersectionIterator &hostIterator )
      : Base( inside, hostIterator )
    {}
  };



  // GeometryGridLevelIntersectionIteratorTraits
  // -------------------------------------------

  template< class Grid >
  struct GeometryGridLevelIntersectionIteratorTraits
  {
    typedef typename remove_const< Grid > :: type :: Traits GridTraits;

    typedef typename GridTraits :: LevelIntersection Intersection;

    typedef typename GridTraits :: HostGrid :: Traits :: LevelIntersectionIterator
    HostIntersectionIterator;
  };



  // GeometryGridLevelIntersectionIterator
  // -------------------------------------

  template< class Grid >
  class GeometryGridLevelIntersectionIterator
    : public GeometryGridIntersectionIterator
      < GeometryGridLevelIntersectionIteratorTraits< Grid > >
  {
    typedef GeometryGridLevelIntersectionIteratorTraits< Grid > Traits;
    typedef GeometryGridIntersectionIterator< Traits > Base;

    typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

  public:
    typedef typename Traits :: Intersection Intersection;

  public:
    template< class Entity >
    GeometryGridLevelIntersectionIterator
      ( const Entity &inside,
      const HostIntersectionIterator &hostIterator )
      : Base( inside, hostIterator )
    {}
  };

}

#endif
