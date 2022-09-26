// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
#define DUNE_PYTHON_GRID_LOCALVIEW_HH

#include <map>

#include <dune/common/visibility.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // LocalViewRegistry
      // -----------------

      template< class LocalView >
      struct DUNE_PRIVATE LocalViewRegistry
      {
        ~LocalViewRegistry()
        {
          for (auto i = binds_.begin(), last = binds_.end(); i != last; ++i)
            i->second = pybind11::object();
          binds_.clear();
        }
        template <class Context>
        void bind ( pybind11::handle localView, pybind11::object context )
        {
          localView.template cast< LocalView & >().bind( context.template cast< const Context & >() );
          find( localView ) = context;
        }

        void unbind ( pybind11::handle localView )
        {
          localView.template cast< LocalView & >().unbind();
          find( localView ) = pybind11::object();
        }

      private:
        pybind11::object &find ( pybind11::handle localView )
        {
          auto result = binds_.emplace( localView.ptr(), pybind11::object() );
          const auto pos = result.first;
          if( result.second )
          {
            pybind11::cpp_function remove_bind( [ this, pos ] ( pybind11::handle weakref ) {
                binds_.erase( pos );
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( localView, remove_bind );
            weakref.release();
          }
          return pos->second;
        }

        std::map< void *, pybind11::object > binds_;
      };



      // localViewRegsitry
      // -----------------

      template< class LocalView >
      inline LocalViewRegistry< LocalView > &localViewRegistry (pybind11::handle self)
      {
        if (!pybind11::hasattr(self, "registry_"))
        {
          auto ptr = new LocalViewRegistry<LocalView>();
          pybind11::cpp_function storage_cleanup(
            [ptr](pybind11::handle weakref) {
              delete ptr;
              weakref.dec_ref();
            });
          (void) pybind11::weakref(self, storage_cleanup).release();
          self.attr("registry_") = (void*)ptr;
        }
        pybind11::handle l = self.attr("registry_");
        return *static_cast< LocalViewRegistry<LocalView>* >(l.cast<void*>());
      }
    } // namespace detail



    // registerLocalView
    // -----------------

    template< class Context, class LocalView, class... options >
    void registerLocalView ( pybind11::class_< LocalView, options... > cls )
    {
      using pybind11::operator""_a;

      cls.def( "bind", [ ] ( pybind11::handle self, pybind11::object context ) {
        detail::localViewRegistry<LocalView>(self).template
                   bind<Context>( self, context ); }, "context"_a );
      cls.def( "unbind", [ ] ( pybind11::handle self ) {
        detail::localViewRegistry<LocalView>(self).unbind( self ); } );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
