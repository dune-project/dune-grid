// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
#define DUNE_GEOGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/geometrygrid/entitypointer.hh>
#include <dune/grid/geometrygrid/intersection.hh>

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
      typedef typename Traits::HostIntersectionIterator HostIntersectionIterator;

      typedef typename Traits::GridTraits GridTraits;

    public:
      typedef typename Traits::Intersection Intersection;
      typedef typename GridTraits::Grid Grid;

      typedef typename GridTraits::template Codim< 0 >::EntityPointer EntityPointer;

    private:
      typedef typename Traits::IntersectionImpl IntersectionImpl;

      typedef typename GridTraits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
      typedef typename GridTraits::template Codim< 0 >::Geometry ElementGeometry;

    public:
      template< class Entity >
      IntersectionIterator ( const Entity &inside,
                             const HostIntersectionIterator &hostIterator )
        : grid_( &inside.grid() ),
          insideGeo_( Grid::getRealImplementation( inside.geometry() ) ),
          hostIterator_( hostIterator ),
          intersection_( 0 )
      {}

      IntersectionIterator  ( const IntersectionIterator &other )
        : grid_( other.grid_ ),
          insideGeo_( Grid::getRealImplementation( other.insideGeo_ ) ),
          hostIterator_( other.hostIterator_ ),
          intersection_( 0 )
      {}

      IntersectionIterator &operator= ( const IntersectionIterator &other )
      {
        grid_ = other.grid_;
        Grid::getRealImplementation( insideGeo_ ) = Grid::getRealImplementation( other.insideGeo_ );
        hostIterator_ = other.hostIterator_;
        update();
        return *this;
      }

      ~IntersectionIterator ()
      {
        if( intersection_ != 0 )
          intersection_->~Intersection();
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
          intersection_ = new( &intersectionMemory_ )Intersection( IntersectionImpl( grid(), insideGeo_, *hostIterator_ ) );
        return *intersection_;
      }

    private:
      const Grid &grid () const
      {
        return *grid_;
      }

      void update ()
      {
        if( intersection_ != 0 )
          intersection_->~Intersection();
        intersection_ = 0;
      }

      const Grid *grid_;
      ElementGeometry insideGeo_;
      HostIntersectionIterator hostIterator_;
      mutable Intersection *intersection_;
      mutable char intersectionMemory_[ sizeof( Intersection ) ];
    };



    // LeafIntersectionIteratorTraits
    // ------------------------------

    template< class Grid >
    struct LeafIntersectionIteratorTraits
    {
      typedef typename remove_const< Grid >::type::Traits GridTraits;

      typedef typename GridTraits::LeafIntersection Intersection;
      typedef LeafIntersection< const Grid > IntersectionImpl;

      typedef typename GridTraits::HostGrid::Traits::LeafIntersectionIterator
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
      typedef typename remove_const< Grid >::type::Traits GridTraits;

      typedef typename GridTraits::LevelIntersection Intersection;
      typedef LevelIntersection< const Grid > IntersectionImpl;

      typedef typename GridTraits::HostGrid::Traits::LevelIntersectionIterator
      HostIntersectionIterator;
    };



    // LevelIntersectionIterator
    // -------------------------

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

#endif // #ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
