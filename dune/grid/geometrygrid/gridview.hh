// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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

    template< class HGV, class CoordFunction, class Allocator >
    class GridView;



    // GridViewTraits
    // --------------

    template< class HGV, class CoordFunction, class Allocator >
    class GridViewTraits
    {
      friend class GridView< HGV, CoordFunction, Allocator >;

      typedef HGV HostGridView;

      typedef typename HostGridView::Grid HostGrid;
      typedef typename HostGridView::Intersection HostIntersection;
      typedef typename HostGridView::IntersectionIterator HostIntersectionIterator;

    public:
      typedef GridView< HostGridView, CoordFunction, Allocator > GridViewImp;

      typedef Dune::GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;

      typedef GeoGrid::IndexSet< const Grid, typename HostGridView::IndexSet > IndexSet;

      typedef Dune::Intersection< const Grid, GeoGrid::Intersection< const Grid, HostIntersection > > Intersection;

      typedef Dune::IntersectionIterator
      < const Grid, GeoGrid::IntersectionIterator< const Grid, HostIntersectionIterator >, GeoGrid::Intersection< const Grid, HostIntersection > >
      IntersectionIterator;

      typedef typename HostGridView::Communication Communication;

      /**
       * \deprecated Use Communication instead! Will be removed after Dune 2.9.
       */
      [[deprecated("Use Communication instead!!")]]
      typedef Communication CollectiveCommunication;

      template< int codim >
      struct Codim
      {
        typedef GeoGrid::Iterator< HostGridView, codim, All_Partition, const Grid > IteratorImp;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImp > Iterator;

        typedef typename Grid::Traits::template Codim< codim >::Entity Entity;

        typedef typename Grid::template Codim< codim >::Geometry Geometry;
        typedef typename Grid::template Codim< codim >::LocalGeometry LocalGeometry;

        template< PartitionIteratorType pit >
        struct Partition
        {
          typedef GeoGrid::Iterator< HostGridView, codim, pit, const Grid > IteratorImp;
          typedef Dune::EntityIterator< codim, const Grid, IteratorImp > Iterator;
        };
      };

      static const bool conforming = HostGridView::conforming;
    };



    // GridView
    // --------

    template< class HGV, class CoordFunction, class Allocator >
    class GridView
    {
      typedef GridView< HGV, CoordFunction, Allocator > This;

    public:
      typedef GridViewTraits< HGV, CoordFunction, Allocator > Traits;

      typedef typename Traits::HostGridView HostGridView;

      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;

      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::Communication Communication;

      /**
       * \deprecated Use Communication instead! Will be removed after Dune 2.9.
       */
      [[deprecated("Use Communication instead!!")]]
      typedef Communication CollectiveCommunication;

      template< int codim >
      struct Codim
        : public Traits::template Codim< codim >
      {};

      static const bool conforming = Traits::conforming;

      GridView ( const Grid &grid, const HostGridView &hostGridView )
        : grid_( &grid ), hostGridView_( hostGridView )
      {}

      GridView ( const This &other )
        : grid_( other.grid_ ), hostGridView_( other.hostGridView_ )
      {}

      GridView ( This &&other )
        : grid_( other.grid_ ), hostGridView_( std::move( other.hostGridView_ ) )
      {}

      This &operator= ( const This &other )
      {
        grid_ = other.grid_;
        hostGridView_ = other.hostGridView_;
        if( indexSet_ )
          indexSet_.reset( hostGridView().indexSet() );
        return *this;
      }

      This &operator= ( This &&other )
      {
        grid_ = other.grid_;
        hostGridView_ = std::move( other.hostGridView_ );
        if( indexSet_ )
          indexSet_.reset( hostGridView().indexSet() );
        return *this;
      }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        indexSet_.reset( hostGridView().indexSet() );
        return indexSet_;
      }

      bool isConforming() const { return hostGridView().isConforming(); }

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
        return begin< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator begin () const
      {
        return Traits::template Codim< codim >::template Partition< pit >::IteratorImp::begin( grid(), hostGridView() );
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        return end< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator end () const
      {
        return Traits::template Codim< codim >::template Partition< pit >::IteratorImp::end( grid(), hostGridView() );
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef GeoGrid::IntersectionIterator< const Grid, typename HostGridView::IntersectionIterator > IntersectionIteratorImpl;
        return IntersectionIteratorImpl( entity, hostGridView().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef GeoGrid::IntersectionIterator< const Grid, typename HostGridView::IntersectionIterator > IntersectionIteratorImpl;
        return IntersectionIteratorImpl( entity, hostGridView().iend( entity.impl().hostEntity() ) );
      }

      const Communication &comm () const
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
      auto communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction ) const
      {
        typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
        typedef GeoGrid::CommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

        WrappedDataHandle wrappedDataHandle( grid(), dataHandle );
        return hostGridView().communicate( wrappedDataHandle, interface, direction );
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
