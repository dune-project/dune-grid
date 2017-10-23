// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_GRID_MAPPER_HH
#define DUNE_PYTHON_GRID_MAPPER_HH

#include <functional>

#include <dune/common/visibility.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>

namespace Dune
{
  namespace Python
  {

    // registerMemberFunctions
    // -----------------------

    template<class GridView, class Mapper, int codim>
    void registerMemberFunctions_(pybind11::class_<Mapper>& cls)
    {
      cls.def("index", &Mapper::template index<typename GridView::template Codim<codim>::Entity>);

      cls.def("contains",
        [] (const Mapper& mapper, const typename GridView::template Codim<codim>::Entity& e) {
          typename Mapper::Index index;
          unsigned res = mapper.contains(e, index);

          return pybind11::make_tuple(res, index);
        });
    }

    template<class GridView, class Mapper, int... codim>
    void registerMemberFunctions(pybind11::class_<Mapper>& cls, std::integer_sequence<int, codim...>)
    {
      std::ignore = std::make_tuple((registerMemberFunctions_<GridView, Mapper, codim>(cls), 0)...);
    }


    // registerMapper
    // --------------

    template<class GridView, class Mapper>
    void registerMapper(pybind11::class_<Mapper> cls)
    {
      // typedef typename Mapper::GridView GridView;
      cls.def( "__len__", [] ( const Mapper &self ) { return self.size(); } );
      cls.def("subIndex", &Mapper::subIndex);
      cls.def("contains",
          [] (Mapper& instance, const typename GridView::template Codim<0>::Entity& e, int i, int cc) {
            typename Mapper::Index index;
            unsigned res = instance.contains(e, i, cc, index);

            return pybind11::make_tuple(res, index);
          });

      cls.def("__call__", [] ( const Mapper &mapper, const typename GridView::template Codim<0>::Entity &element ) {
            // need a cache gt(cdim=0) -> nof indices then we could store directly in retArray
            std::vector<typename Mapper::Index> indices;
            for ( int c=0; c <= GridView::dimension; ++c )
              for ( auto se : range(element.subEntities(c)) )
                for ( auto i : mapper.indices(element, se, c) )
                  indices.push_back(i);
            pybind11::array_t< std::size_t > retArray( indices.size() );
            auto y = retArray.template mutable_unchecked< 1 >();
            std::size_t idx = 0;
            for ( auto i : indices )
              y[idx++] = i;
            return retArray;
          } );

      registerMemberFunctions<GridView, Mapper>(
          cls, std::make_integer_sequence<int, GridView::dimension+1>());
    }

    // registerMultipleCodimMultipleGeomTypeMapper
    // -------------------------------------------

    template<typename GridView>
    auto registerMultipleCodimMultipleGeomTypeMapper(pybind11::handle scope)
    {
      typedef MultipleCodimMultipleGeomTypeMapper<GridView> MCMGMapper;
      auto cls = insertClass<MCMGMapper>(scope, "MultipleCodimMultipleGeomTypeMapper",
          GenerateTypeName("Dune::MultipleCodimMultipleGeomTypeMapper", MetaType<GridView>()),
          IncludeFiles{"dune/grid/common/mcmgmapper.hh","dune/python/grid/mapper.hh"}).first;
      registerMapper<GridView>(cls);
      cls.def( pybind11::init( [] ( const GridView& gridView, pybind11::function contains ) {
            return MCMGMapper(gridView,
                  [contains](Dune::GeometryType gt, int griddim)
                  { return (unsigned)(contains(gt).cast<int>()); }
                ); } ),
          pybind11::keep_alive<1, 2>() );
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_MAPPER_HH
