// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_INDEXSET_HH
#define DUNE_PYTHON_GRID_INDEXSET_HH

#include <dune/common/visibility.hh>
#include <dune/common/version.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace CorePy
  {

    namespace detail
    {

      // registerSubIndex
      // ----------------

      template< class Entity, class IndexSet, class... options >
      inline static std::enable_if_t< Entity::codimension == 0 >
      registerSubIndex ( pybind11::class_< IndexSet, options... > &cls )
      {
        cls.def( "subIndex", [] ( const IndexSet &indexSet, const Entity &entity, int i, int codim ) {
            if( (codim < Entity::codimension) || (codim > Entity::dimension) )
              throw std::invalid_argument( "Invalid codimension: " + std::to_string( codim ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            return static_cast< int >( indexSet.subIndex( entity, i, codim ) );
          } );
      }

      template< class Entity, class Cls >
      inline static void registerSubIndex ( Cls &cls )
      {}



      // registerIndexMethods
      // --------------------

      template< class Entity, class IndexSet, class... options>
      inline static void registerIndexMethods ( pybind11::class_< IndexSet, options... > &cls )
      {
        cls.def( "index", [] ( const IndexSet &indexSet, const Entity &entity ) { return static_cast< int >( indexSet.index( entity ) ); } );
        cls.def( "contains", [] ( const IndexSet &indexSet, const Entity &entity ) { return indexSet.contains( entity ); } );
        registerSubIndex< Entity >( cls );
      }



      // registerGridIndexSet
      // --------------------

      template< class IndexSet, class... Entity >
      inline static void registerGridIndexSet ( pybind11::handle scope, pybind11::class_<IndexSet> cls )
      {
        cls.def( "size", [] ( IndexSet &indexSet, Dune::GeometryType type ) { return indexSet.size( type ); } );
        cls.def( "size", [] ( IndexSet &indexSet, int codim ) { return indexSet.size( codim ); } );

        cls.def( "types", [] ( const IndexSet &indexSet, int codim ) {
            pybind11::list types;
            for( GeometryType type : indexSet.types( codim ) )
              types.append( pybind11::cast( type ) );
            return types;
          } );

        std::ignore = std::make_tuple( (registerIndexMethods< Entity >( cls ), Entity::codimension)... );
      }

    } // namespace detail



    // registerGridIndexSet
    // --------------------

    template< class GridView, class... Entity >
    pybind11::class_< typename GridView::IndexSet > registerGridIndexSet ( pybind11::handle scope )
    {
      auto cls = insertClass< typename GridView::IndexSet >(scope, "IndexSet",
          GenerateTypeName(MetaType<GridView>(),"IndexSet") ).first;
      detail::registerGridIndexSet< typename GridView::IndexSet, Entity... >( scope, cls );
      return cls;
    }



    // registerGridViewIndexSet
    // ------------------------

    template< class GridView, int... codim >
    inline static pybind11::class_< typename GridView::IndexSet > registerGridViewIndexSet ( pybind11::handle scope, std::integer_sequence< int, codim... > )
    {
      return registerGridIndexSet< GridView, typename GridView::template Codim< codim >::Entity... >( scope );
    };

    template< class GridView >
    inline static pybind11::class_< typename GridView::IndexSet > registerGridViewIndexSet ( pybind11::handle scope )
    {
      return registerGridViewIndexSet< GridView >( scope, std::make_integer_sequence< int, GridView::dimension+1 >() );
    };

  } // namespace CorePy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_INDEXSET_HH
