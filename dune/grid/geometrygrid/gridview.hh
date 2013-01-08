// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRIDVIEW_HH
#define DUNE_GEOGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/geometrygrid/datahandle.hh>
#include <dune/grid/geometrygrid/indexsets.hh>
#include <dune/grid/geometrygrid/intersection.hh>
#include <dune/grid/geometrygrid/intersectioniterator.hh>
#include <dune/grid/geometrygrid/iterator.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class HGV, class CoordFunction, class Allocator, PartitionIteratorType pitype >
    class GridView;



    // GridViewTraits
    // --------------

    template< class HGV, class CoordFunction, class Allocator, PartitionIteratorType pitype >
    class GridViewTraits
    {
      friend class GridView< HGV, CoordFunction, Allocator, pitype >;

      typedef HGV HostGridView;

      typedef typename HostGridView::Grid HostGrid;
      typedef typename HostGridView::Intersection HostIntersection;
      typedef typename HostGridView::IntersectionIterator HostIntersectionIterator;

    public:
      typedef GridView< HostGridView, CoordFunction, Allocator, pitype > GridViewImp;

      typedef Dune::GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;

      typedef GeoGrid::IndexSet< const Grid, typename HostGridView::IndexSet > IndexSet;

      typedef Dune::Intersection< const Grid, GeoGrid::Intersection< const Grid, HostIntersection > > Intersection;

      typedef Dune::IntersectionIterator
      < const Grid, GeoGrid::IntersectionIterator< const Grid, HostIntersectionIterator >, GeoGrid::Intersection< const Grid, HostIntersection > >
      IntersectionIterator;

      typedef typename HostGridView::CollectiveCommunication CollectiveCommunication;

      template< int codim >
      struct Codim
      {
        typedef GeoGrid::IteratorTraits< HostGridView, codim, pitype, const Grid > IteratorTraits;
        typedef Dune::EntityIterator< codim, const Grid, GeoGrid::Iterator< IteratorTraits > > Iterator;

        typedef typename Grid::Traits::template Codim< codim >::Entity Entity;
        typedef typename Grid::Traits::template Codim< codim >::EntityPointer EntityPointer;

        typedef typename Grid::template Codim< codim >::Geometry Geometry;
        typedef typename Grid::template Codim< codim >::LocalGeometry LocalGeometry;

        template< PartitionIteratorType pit >
        struct Partition
        {
          typedef GeoGrid::IteratorTraits< HostGridView, codim, pit, const Grid > IteratorTraits;
          typedef Dune::EntityIterator< codim, const Grid, GeoGrid::Iterator< IteratorTraits > > Iterator;
        };
      };

      static const bool conforming = HostGridView::conforming;
    };



    // GridView
    // --------

    template< class HGV, class CoordFunction, class Allocator, PartitionIteratorType pitype >
    class GridView
    {
      typedef GridView< HGV, CoordFunction, Allocator, pitype > This;

    public:
      typedef GridViewTraits< HGV, CoordFunction, Allocator, pitype > Traits;

      typedef typename Traits::HostGridView HostGridView;

      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;

      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< int codim >
      struct Codim
        : public Traits::template Codim< codim >
      {};

      static const bool conforming = Traits::conforming;

      GridView ( const Grid &grid, const HostGridView &hostGridView )
        : grid_( &grid ),
          hostGridView_( hostGridView )
      {}

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        if( !indexSet_ )
          indexSet_ = IndexSet( hostGridView().indexSet() );
        return indexSet_;
      }

      int size ( int codim ) const
      {
        return hostGridView().size( codim );
      }

      int size ( const GeometryType &type ) const
      {
        return hostGridView().size( type );
      }

      template< int codim >
      typename Codim< codim >::Iterator begin () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorTraits IteratorTraits;
        return GeoGrid::Iterator< IteratorTraits >( grid(), hostGridView(), IteratorTraits::begin );
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator begin () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorTraits IteratorTraits;
        return GeoGrid::Iterator< IteratorTraits >( grid(), hostGridView(), IteratorTraits::begin );
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorTraits IteratorTraits;
        return GeoGrid::Iterator< IteratorTraits >( grid(), hostGridView(), IteratorTraits::end );
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator end () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorTraits IteratorTraits;
        return GeoGrid::Iterator< IteratorTraits >( grid(), hostGridView(), IteratorTraits::end );
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef GeoGrid::IntersectionIterator< const Grid, typename HostGridView::IntersectionIterator > IntersectionIteratorImpl;
        return IntersectionIteratorImpl( entity, hostGridView().ibegin( Grid::getRealImplementation( entity ).hostEntity() ) );
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef GeoGrid::IntersectionIterator< const Grid, typename HostGridView::IntersectionIterator > IntersectionIteratorImpl;
        return IntersectionIteratorImpl( entity, hostGridView().iend( Grid::getRealImplementation( entity ).hostEntity() ) );
      }

      const CollectiveCommunication &comm () const
      {
        return hostGridView().comm();
      }

      int overlapSize ( int codim ) const
      {
        return hostGridView().overlapSize( codim );
      }

      int ghostSize ( int codim ) const
      {
        return hostGridView().ghostSize( codim );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction ) const
      {
        typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
        typedef GeoGrid::CommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

        WrappedDataHandle wrappedDataHandle( grid(), dataHandle );
        hostGridView().communicate( wrappedDataHandle, interface, direction );
      }

      const HostGridView &hostGridView () const { return hostGridView_; }

    private:
      const Grid *grid_;
      HostGridView hostGridView_;
      mutable IndexSet indexSet_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GRIDVIEW_HH
