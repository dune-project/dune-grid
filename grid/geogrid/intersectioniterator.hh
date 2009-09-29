// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geogrid/entitypointer.hh>
#include <dune/grid/geogrid/storage.hh>
#include <dune/grid/geogrid/intersection.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid >
    class LeafIntersectionIterator;

    template< class Grid >
    class LevelIntersectionIterator;



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
      typedef GeoGrid :: Storage< IntersectionWrapper > IntersectionStorage;

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
