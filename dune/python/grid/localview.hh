#ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
#define DUNE_PYTHON_GRID_LOCALVIEW_HH

#include <map>

#include <dune/common/visibility.hh>
#include <dune/python/grid/singleton.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // LocalViewRegistry
      // -----------------


      template< class LocalView, class Context >
      struct DUNE_PRIVATE LocalViewRegistry
      {
        ~LocalViewRegistry()
        {
          for (auto i = binds_.begin(), last = binds_.end(); i != last; ++i)
            i->second = pybind11::object();
          binds_.clear();
        }
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

      template< class LocalView, class Context >
      inline LocalViewRegistry< LocalView, Context > &localViewRegistry ()
      {
        static SingletonStorage &singleton = pybind11::cast< SingletonStorage & >( pybind11::module::import( "dune.grid" ).attr( "singleton" ) );
        static LocalViewRegistry< LocalView, Context > &instance = singleton.template instance< LocalViewRegistry<LocalView,Context> >();
        return instance;
      }

    } // namespace detail



    // registerLocalView
    // -----------------

    template< class Context, class LocalView, class... options >
    void registerLocalView ( pybind11::class_< LocalView, options... > cls )
    {
      using pybind11::operator""_a;

      auto &registry = detail::localViewRegistry< LocalView, Context >();
      cls.def( "bind", [ &registry ] ( pybind11::handle self, pybind11::object context ) { registry.bind( self, context ); }, "context"_a );
      cls.def( "unbind", [ &registry ] ( pybind11::handle self ) { registry.unbind( self ); } );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
