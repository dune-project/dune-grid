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

    namespace detail
    {

      // mapperSubIndex
      // --------------

      template< class Mapper, class Entity >
      inline static pybind11::object mapperSubIndex ( const Mapper &mapper, const Entity &entity, int i, int c )
      {
        if( (c < Entity::codimension) || (c > Entity::dimension) )
          throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
        const int size = entity.subEntities( c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        typename Mapper::Index index;
        return (mapper.contains( entity, i, c, index ) ? pybind11::cast( index ) : pybind11::none());
      }



      // registerMapperSubIndex
      // ----------------------

      template< class Entity, class Mapper, class... options >
      inline static std::enable_if_t< Entity::codimension == 0 >
      registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        pybind11::options opts;
        opts.disable_function_signatures();

        cls.def( "subIndex", [] ( const Mapper &self, const Entity &entity, int i, int c ) {
            return detail::mapperSubIndex( self, entity, i, c );
          } );
        cls.def( "subIndex", [] ( const Mapper &self, const Entity &entity, std::tuple< int, int > e ) {
            return detail::mapperSubIndex( self, entity, std::get< 0 >( e ), std::get< 1 >( e ) );
          } );

        cls.def( "subIndices", [] ( const Mapper &self, const Entity &entity, int c ) {
            if( (c < Entity::codimension) || (c > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const int size = entity.subEntities( c );
            pybind11::tuple subIndices( size );
            for( int i = 0; i < size; ++i )
            {
              typename Mapper::Index index;
              subIndices[ i ] = (self.contains( entity, i, c, index ) ? pybind11::cast( index ) : pybind11::none());
            }
            return subIndices;
          }, "entity"_a, "codim"_a );
      }

      template< class Entity, class Mapper, class... options >
      inline static void registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls, PriorityTag< 0 > )
      {}

      template< class Entity, class Mapper, class... options >
      inline static void registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls )
      {
        registerMapperSubIndex< Entity >( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerMapper
    // --------------

    template< class GridView, class Mapper, class... options >
    inline static void registerMapper ( pybind11::class_< Mapper, options... > cls )
    {
      cls.def( "__len__", [] ( const Mapper &self ) { return self.size(); } );
      cls.def_property_readonly( "size", [] ( const Mapper &self ) { return self.size(); } );

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &cls ] ( auto &&codim ) {
          typedef typename GridView::template Codim< codim >::Entity Entity;

          using pybind11::operator""_a;

          cls.def( "index", [] ( const Mapper &self, const Entity &e ) -> pybind11::object {
              typename Mapper::Index index;
              if( self.contains( e, index ) )
                return pybind11::cast( index );
              else
                return pybind11::none();
            }, "entity"_a );
          cls.def( "contains", [] ( const Mapper &self, const Entity &e ) {
              typename Mapper::Index index;
              return self.contains( e, index );
            }, "entity"_a );
          detail::registerMapperSubIndex< Entity >( cls );
        } );

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
