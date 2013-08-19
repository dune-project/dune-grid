// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERING_GRIDVIEW_HH
#define DUNE_GRID_UTILITY_FILTERING_GRIDVIEW_HH

#include <dune/grid/common/entityiterator.hh>

#include "filteringentityiterator.hh"

namespace Dune {

  template<class Filter, class HostView>
  class FilteringGridViewImpl
  {
    Filter filter_;
    HostView host_;

  public:
    FilteringGridViewImpl(const Filter &filter, const HostView &host) :
      filter_(filter), host_(host)
    { }

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const typename HostView::Grid &grid () const
    {
      return host_.grid();
    }

    /** \brief obtain the index set */
    const typename HostView::IndexSet &indexSet () const
    {
      return host_.indexSet();
    }

    // /** \brief obtain number of entities in a given codimension */
    // Not Implemented
    // int size ( int codim ) const
    // {
    //   return host_.size( codim );
    // }

    // /** \brief obtain number of entities with a given geometry type */
    // Not Implemented
    // int size ( const GeometryType &type ) const
    // {
    //   return host_.size( type );
    // }

    /** @brief Return true if the given entity is contained in this grid view
     * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return host_.contains(e) && filter_.contains(e);
    }

    /** \brief obtain begin iterator for this view */
    template<int cd>
    FilteringEntityIteratorImpl<
      Filter,
      typename HostView::template Codim<cd>::Iterator>
    begin() const
    {
      typedef typename HostView::template Codim<cd>::Iterator HI;
      typedef FilteringEntityIteratorImpl<Filter, HI> II;
      return II(filter_, host_.template begin<cd>(),
                host_.template end<cd>());
    }

    /** \brief obtain end iterator for this view */
    template<int cd>
    FilteringEntityIteratorImpl<
      Filter,
      typename HostView::template Codim<cd>::Iterator>
    end() const
    {
      typedef typename HostView::template Codim<cd>::Iterator HI;
      typedef FilteringEntityIteratorImpl<Filter, HI> II;
      return II(filter_, host_.template end<cd>());
    }

    /** \brief obtain begin iterator for this view */
    template<int cd, PartitionIteratorType pit>
    FilteringEntityIteratorImpl<
      Filter,
      typename HostView::template Codim<cd>::template Partition<pit>::
        Iterator>
    begin() const
    {
      typedef typename HostView::template Codim<cd>::
        template Partition<pit>::Iterator HI;
      typedef FilteringEntityIteratorImpl<Filter, HI> II;
      return II(filter_, host_.template begin<cd, pit>(),
                host_.template end<cd, pit>());
    }

    /** \brief obtain end iterator for this view */
    template<int cd, PartitionIteratorType pit>
    FilteringEntityIteratorImpl<
      Filter,
      typename HostView::template Codim<cd>::template Partition<pit>::
        Iterator>
    end() const
    {
      typedef typename HostView::template Codim<cd>::
        template Partition<pit>::Iterator HI;
      typedef FilteringEntityIteratorImpl<Filter, HI> II;
      return II(filter_, host_.template end<cd, pit>());
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    template<class Entity>
    typename HostView::IntersectionIterator
    ibegin(const Entity &entity ) const
    {
      return host_.ibegin(entity);
    }

    /** \brief obtain end intersection iterator with respect to this view */
    template<class Entity>
    typename HostView::IntersectionIterator
    iend ( const Entity &entity ) const
    {
      return host_.iend(entity);
    }

    // /** \brief obtain collective communication object */
    // Not Implemented
    // const CollectiveCommunication &comm () const
    // {
    //   return host_.comm();
    // }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
      return host_.overlapSize(codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
      return host_.ghostSize(codim);
    }

    // /** communicate data on this view */
    // Not Implemented
    // template< class DataHandleImp, class DataType >
    // void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
    //                    InterfaceType iftype,
    //                    CommunicationDirection dir ) const
    // {
    //   host_.communicate(data,iftype,dir);
    // }
  };

  template<class Filter, class HostView>
  struct FilteringGridViewTraits
  {
    typedef FilteringGridViewImpl<Filter, HostView> GridViewImp;

    typedef typename HostView::Grid Grid;
    typedef typename HostView::IndexSet IndexSet;
    typedef typename HostView::Intersection Intersection;
    typedef typename HostView::IntersectionIterator IntersectionIterator;
    typedef typename HostView::CollectiveCommunication
      CollectiveCommunication;

    template<int cd>
    struct Codim
    {
      typedef EntityIterator<
        cd, typename HostView::Grid,
        FilteringEntityIteratorImpl<
          Filter,
          typename HostView::template Codim<cd>::Iterator> > Iterator;
      typedef typename HostView::template Codim<cd>::EntityPointer
        EntityPointer;
      typedef typename HostView::template Codim<cd>::Entity Entity;
      typedef typename HostView::template Codim<cd>::Geometry Geometry;
      typedef typename HostView::template Codim<cd>::LocalGeometry
        LocalGeometry;

      template< PartitionIteratorType pit >
      struct Partition
      {
      typedef EntityIterator<
        cd, typename HostView::Grid,
        FilteringEntityIteratorImpl<
          Filter,
          typename HostView::template Codim<cd>::
            template Partition<pit>::Iterator> > Iterator;
      };
    };

    enum { conforming = HostView::conforming };
    typedef typename HostView::ctype ctype;
    enum { dimension = HostView::dimension };
    enum { dimensionworld = HostView::dimensionworld };

  };

  template<class Filter, class HostView>
  class FilteringGridView :
    public GridView<FilteringGridViewTraits<Filter, HostView> >
  {
    typedef FilteringGridViewTraits<Filter, HostView> T;
    typedef typename T::GridViewImp Impl;
    typedef GridView<T> Facade;
  public:
    FilteringGridView(const Filter &filter, const HostView host) :
      Facade(Impl(filter, host))
    { }

    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return this->impl().contains(e);
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_FILTERING_GRIDVIEW_HH
