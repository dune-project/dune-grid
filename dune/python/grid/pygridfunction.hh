// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_PY_PYFUNCTION_HH
#define DUNE_GRID_PY_PYFUNCTION_HH

#include <string>
#include <utility>

#include <dune/python/pybind11/pybind11.h>
#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Python
  {

    // PyGridFunction
    // --------------

    template< class GridFunction >
    class DUNE_PRIVATE PyGridFunction
    {
    public:
      PyGridFunction ( const GridFunction &impl, pybind11::object pyObj )
        : pyObj_( std::move( pyObj ) ), lf_(impl)
      {}

      PyGridFunction ( const GridFunction &impl )
        : pyObj_( pybind11::reinterpret_borrow<pybind11::object>(
                  pybind11::detail::get_object_handle( &impl, pybind11::detail::get_type_info( typeid( GridFunction ) ) )
                ) ),
          lf_(impl)
      {}

      template< class Point >
      auto evaluate ( const Point &x ) const { return lf_( x ); }
      template< class Point >
      auto operator() ( const Point &x ) const { return lf_( x ); }
      template <class Entity>
      void bind(const Entity &entity) { lf_.bind(entity); }
      void unbind() { lf_.unbind(); }

      // friend PyGridFunction<GridFunction> localFunction ( const PyGridFunction &gf ) { return PyGridFunction( gf ); }

    protected:
      pybind11::object pyObj_;
      typename GridFunction::LocalFunction lf_;
    };

    // pyGridFunction
    // --------------

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction );
    }

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction, pybind11::object pyObj ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction, std::move( pyObj ) );
    }

  } // namespace Python

} // namespace Dune

#endif // DUNE_GRID_PY_PYFUNCTION_HH
