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
        : hostIterator_( hostIterator ),
          intersection_( IntersectionImpl( inside.geometry() ) )
      {}

      IntersectionIterator ( const IntersectionIterator &other )
        : hostIterator_( other.hostIterator_ ),
          intersection_( IntersectionImpl( Grid::getRealImplementation( other.intersection_ ) ) )
      {}

      IntersectionIterator &operator= ( const IntersectionIterator &other )
      {
        hostIterator_ = other.hostIterator_;
        Grid::getRealImplementation( intersection_ ) = Grid::getRealImplementation( other.intersection_ );
        return *this;
      }

      bool equals ( const IntersectionIterator &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }

      void increment ()
      {
        ++hostIterator_;
        intersectionImpl().invalidate();
      }

      const Intersection &dereference () const
      {
        if( !intersectionImpl() )
          intersectionImpl().initialize( *hostIterator_ );
        return intersection_;
      }

    private:
      IntersectionImpl &intersectionImpl () const
      {
        return Grid::getRealImplementation( intersection_ );
      }

      HostIntersectionIterator hostIterator_;
      mutable Intersection intersection_;
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

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_INTERSECTIONITERATOR_HH
