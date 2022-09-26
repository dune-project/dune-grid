// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_FUNCTION_HH
#define DUNE_PYTHON_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/visibility.hh>

#include <dune/python/common/dimrange.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/common/vector.hh>
#include <dune/python/common/fvector.hh>
#include <dune/python/grid/simplegridfunction.hh>
#include <dune/python/grid/localview.hh>
#include <dune/python/grid/entity.hh>
#include <dune/python/grid/numpy.hh>
#include <dune/python/grid/object.hh>
#include <dune/python/grid/vtk.hh>

#if HAVE_DUNE_VTK
#include <dune/vtk/function.hh>
#include <dune/python/grid/pygridfunction.hh>
#endif

#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // GridFunctionTraits
    // ------------------

    template< class GridFunction >
    struct GridFunctionTraits
      : public GridObjectTraits< GridFunction >
    {
      typedef typename GridObjectTraits< GridFunction >::LocalCoordinate LocalCoordinate;

      typedef std::decay_t< decltype( localFunction( std::declval< const GridFunction & >() ) ) > LocalFunction;
      typedef std::decay_t< decltype( std::declval< LocalFunction & >()( std::declval< const LocalCoordinate & >() ) ) > Range;

      typedef typename GridFunction::GridView GridView;
    };



    namespace detail
    {

      template< class LocalCoordinate, class LocalFunction, class X >
      inline static auto callLocalFunction ( LocalFunction &&f, const X &x, PriorityTag< 2 > )
        -> decltype( f( x ) )
      {
        return f( x );
      }

      template< class LocalCoordinate, class LocalFunction >
      inline static pybind11::object callLocalFunction ( LocalFunction &&f, pybind11::array_t< typename FieldTraits< LocalCoordinate >::field_type > x, PriorityTag< 1 > )
      {
        return vectorize( [ &f ] ( const LocalCoordinate &x ) { return f( x ); }, x );
      }

      template< class LocalCoordinate, class LocalFunction, class X >
      inline static auto callLocalFunction ( LocalFunction &&f, const X &x, PriorityTag<0> )
        -> std::enable_if_t< !std::is_const< std::remove_reference_t< LocalFunction > >::value, pybind11::object >
      {
        return callLocalFunction< LocalCoordinate >( std::forward< LocalFunction >( f ), x, PriorityTag< 42 >() );
      }

      template< class GridFunction, class... options >
      void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
      {
        using pybind11::operator""_a;
        typedef typename GridFunctionTraits< GridFunction >::Range Range;
        cls.def_property_readonly( "grid", [] ( const GridFunction &self ) -> pybind11::handle
                    { return pybind11::cast(Dune::Python::gridView( self )); } );
        cls.def_property_readonly( "dimRange", [] ( pybind11::object self ) { return pybind11::int_( DimRange< Range >::value ); } );
        cls.def( "addToVTKWriter", &addToVTKWriter< GridFunction >, pybind11::keep_alive< 3, 1 >(), "name"_a, "writer"_a, "dataType"_a );

        cls.def( "cellData", [] ( const GridFunction &self, int level ) { return cellData( self, level ); }, "level"_a = 0 );
        cls.def( "pointData", [] ( const GridFunction &self, int level ) { return pointData( self, level ); }, "level"_a = 0 );
        cls.def( "polygonData", [] ( const GridFunction &self ) { return polygonData( self ); },
          R"doc(
            Store the grid with piecewise constant data in numpy arrays.

            Returns: pair with coordinate array storing the vertex coordinate of each polygon
                     in the grid and an array with a range type for each polygon.
          )doc" );
      }
    } // namespace detail

    // registerGridFunction
    // --------------------

    template< class GridFunction, class... options >
    void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
    {
      using pybind11::operator""_a;

      typedef typename GridFunctionTraits< GridFunction >::Element Element;
      typedef typename GridFunctionTraits< GridFunction >::LocalCoordinate LocalCoordinate;
      typedef typename GridFunctionTraits< GridFunction >::LocalFunction LocalFunction;
      typedef typename GridFunctionTraits< GridFunction >::Range Range;

      typedef pybind11::array_t< typename FieldTraits< LocalCoordinate >::field_type > Array;

      // TODO subclassing from a non registered traits class not covered by TypeRegistry
      pybind11::class_< LocalFunction > clsLocalFunction( cls, "LocalFunction", pybind11::dynamic_attr() );
      registerLocalView< Element >( clsLocalFunction );
      clsLocalFunction.def( "__call__", [] ( LocalFunction &self, const LocalCoordinate &x ) {
          return detail::callLocalFunction< LocalCoordinate >( self, x, PriorityTag<2>() );
        }, "x"_a );
      clsLocalFunction.def( "__call__", [] ( LocalFunction &self, Array x ) {
          return detail::callLocalFunction< LocalCoordinate >( self, x, PriorityTag<2>() );
        }, "x"_a );
      clsLocalFunction.def_property_readonly( "dimRange", [] ( pybind11::object self ) { return pybind11::int_( DimRange< Range >::value ); } );

      cls.def( "localFunction", [] ( const GridFunction &self ) { return localFunction( self ); }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "localFunction", [] ( const GridFunction &self, const Element &element )
          { auto lf = localFunction(self); lf.bind(element); return lf; },
          pybind11::keep_alive< 0, 1 >(),
          pybind11::keep_alive< 0, 2 >() );

      cls.def( "__call__", [] ( const GridFunction &self, const Element &element, LocalCoordinate &x ) {
          auto lf = localFunction(self);
          lf.bind(element);
          auto y = detail::callLocalFunction< LocalCoordinate >( lf, x, PriorityTag<2>() );
          lf.unbind();
          return y;
        }, "element"_a, "x"_a );
      cls.def( "__call__", [] ( const GridFunction &self, const Element &element, Array x ) {
          auto lf = localFunction(self);
          lf.bind(element);
          auto y = detail::callLocalFunction< LocalCoordinate >( lf, x, PriorityTag<2>() );
          lf.unbind();
          return y;
        }, "element"_a, "x"_a );
      detail::registerGridFunction(scope,cls);
#if HAVE_DUNE_VTK
      typedef typename GridFunctionTraits< GridFunction >::GridView GridView;
      using VirtualizedGF = Dune::Vtk::Function<GridView>;
      // register the Function class if not already available
      auto vgfClass = Python::insertClass<VirtualizedGF>(scope,"VtkFunction",
          Python::GenerateTypeName("Dune::Vtk::Function",MetaType<GridView>()),
          Python::IncludeFiles{"dune/vtk/function.hh"});
      vgfClass.first.def( pybind11::init( [] ( GridFunction &gf ) {
          // TODO: perhpas grid functions should just have a name attribute in general
          return new VirtualizedGF( pyGridFunction(gf), "tmp" );
        } ) );
      pybind11::implicitly_convertible<GridFunction,VirtualizedGF>();
#endif
    }

    template <class GridView, int dimR>
    struct stdFunction
    {
      static const unsigned int dimRange = (dimR ==0 ? 1 : dimR);
      typedef typename GridView::template Codim< 0 >::Entity Entity;
      typedef typename Entity::Geometry::LocalCoordinate Coordinate;
      typedef typename std::conditional< dimR == 0, double, Dune::FieldVector< double, dimRange > >::type Value;
      typedef std::function<Value(const Entity&,const Coordinate&)> type;
    };
    template <class GridView,int dimR,class Evaluate>
    struct EvaluateType
    {
      typedef typename GridView::template Codim< 0 >::Entity Entity;
      typedef typename Entity::Geometry::LocalCoordinate Coordinate;
      typedef typename std::conditional< dimR == 0, double, Dune::FieldVector< double, dimR > >::type Value;
      static std::string name()
      { std::string entity = findInTypeRegistry<Entity>().first->second.name;
        std::string coord = findInTypeRegistry<Coordinate>().first->second.name;
        std::string value;
        if (dimR==0) value = "double";
        else
        {
          auto found = findInTypeRegistry<Value>();
          assert(!found.second);
          value = found.first->second.name;
        }
        return "std::function<"+value+"(const "+entity+"&,const "+coord+"&)>";
      }
    };
    template <class GridView,int dimR>
    struct EvaluateType<GridView,dimR,pybind11::function>
    {
      static std::string name() { return "pybind11::function"; }
    };

    namespace detail
    {

      // PyGridFunctionEvaluator
      // -----------------------
      template <class GridView, int dimR, class Evaluate>
      struct DUNE_PRIVATE PyGridFunctionEvaluator
      {};

      template <class GridView, int dimR>
      struct DUNE_PRIVATE PyGridFunctionEvaluator<GridView,dimR,pybind11::function>
      {
        static const unsigned int dimRange = (dimR ==0 ? 1 : dimR);

        typedef typename GridView::template Codim< 0 >::Entity Entity;
        typedef typename Entity::Geometry::LocalCoordinate Coordinate;

        typedef typename std::conditional< dimR == 0, double, Dune::FieldVector< double, dimRange > >::type Value;

        explicit PyGridFunctionEvaluator ( pybind11::function evaluate ) : evaluate_( evaluate ) {}

        Value operator() ( const Entity &entity, const Coordinate &x ) const
        {
          pybind11::gil_scoped_acquire acq;
          return pybind11::cast< Value >( evaluate_( entity, x ) );
        }

        pybind11::array_t< double > operator() ( const Entity &entity, const pybind11::array_t<double> x ) const
        {
          pybind11::gil_scoped_acquire acq;
          return pybind11::cast< pybind11::array_t< double > >( evaluate_( entity, x ) );
        }

      private:
        pybind11::function evaluate_;
      };
      template <class GridView, int dimR>
      struct DUNE_PRIVATE PyGridFunctionEvaluator<GridView,dimR,
                          typename stdFunction<GridView,dimR>::type >
      {
        static const unsigned int dimRange = (dimR ==0 ? 1 : dimR);

        typedef typename GridView::template Codim< 0 >::Entity Entity;
        typedef typename Entity::Geometry::LocalCoordinate Coordinate;

        typedef typename std::conditional< dimR == 0, double, Dune::FieldVector< double, dimRange > >::type Value;

        typedef typename stdFunction<GridView,dimR>::type Evaluate;
        explicit PyGridFunctionEvaluator ( Evaluate evaluate ) : evaluate_( evaluate ) {}

        Value operator() ( const Entity &entity, const Coordinate &x ) const
        {
          return evaluate_( entity, x );
        }
      private:
        Evaluate evaluate_;
      };



      // registerPyGridFunction
      // ----------------------

      template< class GridView, class Evaluate, unsigned int dimRange >
      auto registerPyGridFunction ( pybind11::handle scope, const std::string &name, bool scalar, std::integral_constant< unsigned int, dimRange > )
      {
        using pybind11::operator""_a;

        typedef typename GridView::template Codim<0>::Entity Entity;
        typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;
        typedef PyGridFunctionEvaluator<GridView,dimRange,Evaluate> Evaluator;
        typedef SimpleGridFunction< GridView, Evaluator > GridFunction;
        if (dimRange>0)
          Dune::Python::registerFieldVector<double,dimRange>(scope);
        addToTypeRegistry<Evaluator>(GenerateTypeName("Dune::Python::detail::PyGridFunctionEvaluator",
                                            MetaType<GridView>(),dimRange,
                                            EvaluateType<GridView,dimRange,Evaluate>::name()
                                            ),
                         IncludeFiles{"dune/python/grid/function.hh"});

        std::string clsName = name + std::to_string( dimRange );
        auto gf = insertClass< GridFunction >( scope, clsName,
                  pybind11::dynamic_attr(),
                  GenerateTypeName("Dune::Python::SimpleGridFunction",
                                      MetaType<GridView>(), Dune::MetaType<Evaluator>()),
            IncludeFiles{"dune/python/grid/function.hh"});
        gf.first.def(pybind11::init([](GridView &gridView, Evaluate callable) {
              return new GridFunction( gridView,
                         PyGridFunctionEvaluator<GridView,dimRange,Evaluate>(callable) );
              }), "gridView"_a, "callable"_a, pybind11::keep_alive<1,2>() );

        if (gf.second)
        {
          Dune::Python::registerGridFunction( scope, gf.first );
          gf.first.def_property_readonly( "scalar", [scalar] ( pybind11::object self ) { return scalar; } );
          if constexpr (dimRange>0)
          {
            typedef typename Dune::Python::stdFunction<GridView,0>::type Evaluate0;
            detail::registerPyGridFunction< GridView, Evaluate0, 0 >
                    ( scope, name, true, std::integral_constant< unsigned int, 0 >() );
            gf.first.def( "__getitem__", [] ( const GridFunction &self, std::size_t c ) {
              Evaluate0 eval0 = [&self,c](const Entity &e, const LocalCoordinate &x) -> double
              {
                auto lf = localFunction(self);
                lf.bind(e);
                auto y = detail::callLocalFunction< LocalCoordinate >( lf, x, PriorityTag<2>() );
                lf.unbind();
                return y[0];
              };
              auto gridFunction = simpleGridFunction( gridView(self),
                  detail::PyGridFunctionEvaluator< GridView, 0, Evaluate0 >( std::move( eval0 ) )
                );
              return gridFunction;
              // return pybind11::cast( std::move( gridFunction ) );
            }, pybind11::keep_alive< 0, 1 >() );
          }
        }
        return gf;
      }

    } // namespace detail

    template< class GridView, class Evaluate, int dimRange >
    auto registerGridFunction ( pybind11::handle scope, std::string name, bool scalar )
    {
      detail::registerPyGridFunction< GridView, Evaluate, dimRange >( scope, name, scalar, std::integral_constant< unsigned int, dimRange >() );
    }
    template <class Value> struct FunctionRange
    {
      static constexpr int value()
      { if constexpr (std::is_convertible_v<Value,double>) return 0; else return Value::dimension; }
    };
    template< class GridView, class Eval >
    auto registerGridFunction ( pybind11::handle scope, pybind11::object gp, std::string name, Eval eval )
    {
      typedef typename GridView::template Codim<0>::Entity Entity;
      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;
      typedef decltype(eval(std::declval<const Entity&>(),std::declval<const LocalCoordinate&>())) Value;
      static constexpr int dimRange = FunctionRange<Value>::value();
      typedef typename Dune::Python::stdFunction<GridView,dimRange>::type Evaluate;
      registerGridFunction< GridView, Evaluate, dimRange >( scope, name, dimRange==0 );

      Evaluate evaluate(eval);
      const GridView &gridView = gp.cast< const GridView & >();
      auto gridFunction = simpleGridFunction( gridView, detail::PyGridFunctionEvaluator< GridView, dimRange, Evaluate >( std::move( evaluate ) ) );
      return pybind11::cast( std::move( gridFunction ), pybind11::return_value_policy::move, gp );
    }
  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_FUNCTION_HH
