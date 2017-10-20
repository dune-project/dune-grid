// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_GRID_RANGE_HH
#define DUNE_PYTHON_GRID_RANGE_HH

#include <string>
#include <utility>

#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // PyIterator
    // ----------

    template< class Iterator, class Entity >
    struct PyIterator
    {
      PyIterator ( Iterator begin, Iterator end )
        : it_( std::move( begin ) ), end_( std::move( end ) )
      {}

      Entity next ()
      {
        if( it_ == end_ )
          throw pybind11::stop_iteration();

        Entity entity = *it_;
        ++it_;
        return entity;
      }

    private:
      Iterator it_, end_;
    };

    template<typename GridView, int codim>
    using PyGridViewIterator =
      PyIterator<
        typename GridView::template Codim<codim>::Iterator,
        typename GridView::template Codim<codim>::Entity >;

    template<typename GridView, int codim, int partition>
    using PyGridViewParIterator =
      PyIterator<
        typename GridView::template Codim<codim>::
          template Partition< static_cast<Dune::PartitionIteratorType>(partition) >::Iterator,
        typename GridView::template Codim<codim>::Entity >;

    template<typename GridView>
    using PyIntersectionIterator =
      PyIterator<
        typename GridView::IntersectionIterator,
        typename GridView::Intersection >;

    // registerPyIterator
    // ------------------

    template<typename Iterator>
    void registerPyIterator(pybind11::handle scope, pybind11::class_<Iterator> cls)
    {
      cls.def("__iter__", [] (Iterator& it) -> Iterator& { return it; });
      cls.def("__next__", &Iterator::next);
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_RANGE_HH
