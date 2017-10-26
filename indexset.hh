// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_INDEXSET_HH
#define DUNE_PYTHON_GRID_INDEXSET_HH

#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typeutilities.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // indexSetSubIndex
      // ----------------

      template< class IndexSet, class Entity >
      inline static int indexSetSubIndex ( const IndexSet &indexSet, const Entity &entity, int i, int c )
      {
        if( (c < Entity::codimension) || (c > Entity::dimension) )
          throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
        const int size = entity.subEntities( c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        return static_cast< int >( indexSet.subIndex( entity, i, c ) );
      }



      // registerSubIndex
      // ----------------

      template< class Entity, class IndexSet, class... options >
      inline static std::enable_if_t< Entity::codimension == 0 >
      registerSubIndex ( pybind11::class_< IndexSet, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        pybind11::options opts;
        opts.disable_function_signatures();

        cls.def( "subIndex", [] ( const IndexSet &self, const Entity &entity, int i, int c ) {
            return detail::indexSetSubIndex( self, entity, i, c );
          } );
        cls.def( "subIndex", [] ( const IndexSet &self, const Entity &entity, std::tuple< int, int > e ) {
            return detail::indexSetSubIndex( self, entity, std::get< 0 >( e ), std::get< 1 >( e ) );
          } );

        cls.def( "subIndices", [] ( const IndexSet &self, const Entity &entity, int c ) {
            if( (c < Entity::codimension) || (c > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const int size = entity.subEntities( c );
            pybind11::tuple subIndices( size );
            for( int i = 0; i < size; ++i )
              subIndices[ i ] = pybind11::cast( self.subIndex( entity, i, c ) );
            return subIndices;
          }, "entity"_a, "codim"_a );
      }

      template< class Entity, class IndexSet, class... options >
      inline static void registerSubIndex ( pybind11::class_< IndexSet, options... > cls, PriorityTag< 0 > )
      {}

      template< class Entity, class IndexSet, class... options >
      inline static void registerSubIndex ( pybind11::class_< IndexSet, options... > cls )
      {
        registerSubIndex< Entity >( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerGridViewIndexSet
    // ------------------------

    template< class GridView, class IndexSet, class... options >
    inline static void registerGridViewIndexSet ( pybind11::handle scope, pybind11::class_< IndexSet, options... > cls )
    {
      using pybind11::operator""_a;

      pybind11::options opts;
      opts.disable_function_signatures();

      cls.def( "size", [] ( IndexSet &self, Dune::GeometryType type ) { return self.size( type ); } );
      cls.def( "size", [] ( IndexSet &self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])" );
          return self.size( codim );
        } );

      cls.def( "types", [] ( const IndexSet &self, int codim ) {
          pybind11::list types;
          for( GeometryType type : self.types( codim ) )
            types.append( pybind11::cast( type ) );
          return types;
        }, "codim"_a );

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &cls ] ( auto &&codim ) {
          typedef typename GridView::template Codim< codim >::Entity Entity;

          using pybind11::operator""_a;

          cls.def( "index", [] ( const IndexSet &self, const Entity &entity ) { return static_cast< int >( self.index( entity ) ); }, "entity"_a );
          cls.def( "contains", [] ( const IndexSet &self, const Entity &entity ) { return self.contains( entity ); }, "entity"_a );
          detail::registerSubIndex< Entity >( cls );
        } );
    }



    // registerGridViewIndexSet
    // ------------------------

    template< class GridView >
    inline static pybind11::class_< typename GridView::IndexSet > registerGridViewIndexSet ( pybind11::handle scope )
    {
      typedef typename GridView::IndexSet IndexSet;

      auto cls = insertClass< IndexSet >( scope, "IndexSet", GenerateTypeName( MetaType< GridView >(), "IndexSet" ) ).first;
      registerGridViewIndexSet< GridView >( scope, cls );
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_INDEXSET_HH
