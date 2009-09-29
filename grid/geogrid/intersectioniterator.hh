// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geogrid/entitypointer.hh>
#include <dune/grid/geogrid/cornerstorage.hh>
#include <dune/grid/geogrid/storage.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid, class HostIntersection >
    class Intersection;

    template< class Grid >
    class LeafIntersection;

    template< class Grid >
    class LevelIntersection;

    template< class Intersection >
    class IntersectionWrapper;

    template< class Grid >
    class LeafIntersectionIterator;

    template< class Grid >
    class LevelIntersectionIterator;



    // GeometryGridIntersectionIterator
    // --------------------------------

    template< class HostGrid, class CoordFunction, class HostIntersection >
    class Intersection< const GeometryGrid< HostGrid, CoordFunction >, HostIntersection >
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
      typedef GeoGrid :: CoordVector< dimension-1, const Grid, false > CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

      const EntityPointer *inside_;
      const HostIntersection *hostIntersection_;
      mutable MakeableGeometry geo_;

    public:
      Intersection ()
        : geo_( GeometryImpl() )
      {}

      Intersection ( const Intersection &other )
        : inside_( other.inside_ ),
          hostIntersection_( other.hostIntersection_ ),
          geo_( GeometryImpl() )
      {}

      const EntityPointer &inside () const
      {
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


    // LeafIntersection
    // ----------------

    template< class HostGrid, class CoordFunction >
    class LeafIntersection< const GeometryGrid< HostGrid, CoordFunction > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction >,
            typename HostGrid :: Traits :: LeafIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;
      typedef typename HostGrid :: Traits :: LeafIntersection HostIntersection;

      template< class > friend class IntersectionWrapper;
    };



    // LevelIntersection
    // -----------------

    template< class HostGrid, class CoordFunction >
    class LevelIntersection< const GeometryGrid< HostGrid, CoordFunction > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction >,
            typename HostGrid :: Traits :: LevelIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;
      typedef typename HostGrid :: Traits :: LevelIntersection HostIntersection;

      template< class > friend class IntersectionWrapper;
    };



    // GeomegryGridIntersectionWrapper
    // -------------------------------

    template< class Intersection >
    class IntersectionWrapper
      : public Intersection
    {
      typedef Intersection Base;

    protected:
      using Base :: getRealImp;

    public:
      typedef typename Intersection :: ImplementationType Implementation;

      typedef typename Implementation :: EntityPointer EntityPointer;
      typedef typename Implementation :: HostIntersection HostIntersection;

      IntersectionWrapper ()
        : Base( Implementation() )
      {}

      void initialize( const EntityPointer &inside,
                       const HostIntersection &hostIntersection )
      {
        getRealImp().initialize( inside, hostIntersection );
      }
    };



    // IntersectionIterator
    // --------------------

    template< class Traits >
    class IntersectionIterator
    {
      typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

    public:
      typedef typename Traits :: Intersection Intersection;
      typedef typename Traits :: GridTraits :: Grid Grid;

      typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;

    private:
      typedef GeoGrid :: IntersectionWrapper< Intersection > IntersectionWrapper;
      typedef GeometryGridStorage< IntersectionWrapper > IntersectionStorage;

      EntityPointer inside_;
      HostIntersectionIterator hostIterator_;
      mutable IntersectionWrapper *intersection_;

    public:
      template< class Entity >
      IntersectionIterator ( const Entity &inside,
                             const HostIntersectionIterator &hostIterator )
        : inside_( inside.template entity< 0 >( 0 ) ),
          hostIterator_( hostIterator ),
          intersection_( 0 )
      {}

      IntersectionIterator  ( const IntersectionIterator &other )
        : inside_( other.inside_ ),
          hostIterator_( other.hostIterator_ ),
          intersection_( 0 )
      {}

      IntersectionIterator &operator= ( const IntersectionIterator &other )
      {
        inside_ = other.inside_;
        hostIterator_ = other.hostIterator_;
        update();
        return *this;
      }

      ~IntersectionIterator ()
      {
        IntersectionStorage :: free( intersection_ );
      }

      bool equals ( const IntersectionIterator &other ) const
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



    // LeafIntersectionIteratorTraits
    // ------------------------------

    template< class Grid >
    struct LeafIntersectionIteratorTraits
    {
      typedef typename remove_const< Grid > :: type :: Traits GridTraits;

      typedef typename GridTraits :: LeafIntersection Intersection;

      typedef typename GridTraits :: HostGrid :: Traits :: LeafIntersectionIterator
      HostIntersectionIterator;
    };



    // LeafIntersectionIterator
    // ------------------------

    template< class Grid >
    class LeafIntersectionIterator
      : public IntersectionIterator< LeafIntersectionIteratorTraits< Grid > >
    {
      typedef LeafIntersectionIteratorTraits< Grid > Traits;
      typedef IntersectionIterator< Traits > Base;

      typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

    public:
      typedef typename Traits :: Intersection Intersection;

    public:
      template< class Entity >
      LeafIntersectionIterator ( const Entity &inside,
                                 const HostIntersectionIterator &hostIterator )
        : Base( inside, hostIterator )
      {}
    };



    // LevelIntersectionIteratorTraits
    // -------------------------------

    template< class Grid >
    struct LevelIntersectionIteratorTraits
    {
      typedef typename remove_const< Grid > :: type :: Traits GridTraits;

      typedef typename GridTraits :: LevelIntersection Intersection;

      typedef typename GridTraits :: HostGrid :: Traits :: LevelIntersectionIterator
      HostIntersectionIterator;
    };



    // LevelIntersectionIterator
    // -------------------------------------

    template< class Grid >
    class LevelIntersectionIterator
      : public IntersectionIterator< LevelIntersectionIteratorTraits< Grid > >
    {
      typedef LevelIntersectionIteratorTraits< Grid > Traits;
      typedef IntersectionIterator< Traits > Base;

      typedef typename Traits :: HostIntersectionIterator HostIntersectionIterator;

    public:
      typedef typename Traits :: Intersection Intersection;

    public:
      template< class Entity >
      LevelIntersectionIterator ( const Entity &inside,
                                  const HostIntersectionIterator &hostIterator )
        : Base( inside, hostIterator )
      {}
    };

  }

}

#endif
