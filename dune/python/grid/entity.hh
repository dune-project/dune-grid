// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_ENTITY_HH
#define DUNE_PYTHON_GRID_ENTITY_HH

#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/grid/geometry.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/operators.h>

namespace Dune
{

  namespace Python
  {

    // PyHierarchicIterator
    // --------------------

    template<class Entity>
    struct PyHierarchicIterator
    {
      typedef typename Entity::HierarchicIterator HierarchicIterator;

      PyHierarchicIterator(const Entity& e, int maxLevel)
        : it_(e.hbegin(maxLevel)), end_(e.hend(maxLevel))
      {}

      auto next()
      {
        if (it_ == end_)
          throw pybind11::stop_iteration();

        return *it_++;
      }

    private:
      HierarchicIterator it_;
      HierarchicIterator end_;
    };


    // registerPyHierarchicIterator
    // ----------------------------

    template<class Entity>
    auto registerPyHierarchicIterator(pybind11::handle scope)
    {
      typedef PyHierarchicIterator<Entity> Iterator;

      auto cls = insertClass<Iterator>(scope, "PyHierarchicIterator",
          GenerateTypeName("PyHierarchicIterator",MetaType<Entity>()),
          IncludeFiles{"dune/python/entity.hh"}).first;
       cls.def("__iter__", [] (Iterator& it) -> Iterator& { return it; });
       cls.def("__next__", &Iterator::next);
    }



    namespace detail
    {

      // makeSubEntity
      // -------------

      template< class Entity, int codim >
      inline static auto makeSubEntity ( const Entity &entity, int i )
        -> pybind11::object
      {
        const int size = entity.subEntities( codim );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        return pybind11::cast( entity.template subEntity< codim >( i ) );
      }



      // makeSubEntities
      // ---------------

      template< class Entity, int codim >
      inline static auto makeSubEntities ( const Entity &entity )
        -> pybind11::tuple
      {
        const int size = entity.subEntities( codim );
        pybind11::tuple subEntities( size );
        for( int i = 0; i < size; ++i )
          subEntities[ i ] = pybind11::cast( entity.template subEntity< codim >( i ) );
        return subEntities;
      }



      // registerExtendedEntityInterface
      // -------------------------------

      template< class Entity, class... options >
      inline static auto registerExtendedEntityInterface ( pybind11::class_< Entity, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< Entity::codimension == 0 >
      {
        const int dimension = Entity::dimension;
        const int mydimension = Entity::mydimension;
        const int codimension = Entity::codimension;

        static_assert( mydimension + codimension == dimension, "Entity exports incompatible dimensions" );

        using pybind11::operator""_a;

        cls.def_property_readonly( "father", [] ( const Entity &self ) {
            return (self.hasFather() ? pybind11::cast( self.father() ) : pybind11::none());
          } );

        cls.def_property_readonly( "geometryInFather", [] ( const Entity &self ) {
            return (self.hasFather() ? pybind11::cast( self.geometryInFather() ) : pybind11::none());
          } );

        cls.def_property_readonly( "isLeaf", [] ( const Entity &self ) { return self.isLeaf(); } );
        cls.def_property_readonly( "isRegular", [] ( const Entity &self ) { return self.isRegular(); } );
        cls.def_property_readonly( "isNew", [] ( const Entity &self ) { return self.isNew(); } );
        cls.def_property_readonly( "mightVanish", [] ( const Entity &self ) { return self.mightVanish(); } );
        cls.def( "hasBoundaryIntersections", [] ( const Entity &self ) { return self.hasBoundaryIntersections(); } );

        registerPyHierarchicIterator< Entity >( cls );
        cls.def( "descendants", [] ( const Entity &self, int maxLevel ) {
            return PyHierarchicIterator< Entity >( self, maxLevel );
          }, pybind11::keep_alive< 0, 1 >() );

        std::array< pybind11::object (*) ( const Entity &, int ), dimension+1 > makeSubEntity;
        std::array< pybind11::tuple (*) ( const Entity & ), dimension+1 > makeSubEntities;
        Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ &makeSubEntity, &makeSubEntities ] ( auto codim ) {
            makeSubEntity[ codim ] = detail::makeSubEntity< Entity, codim >;
            makeSubEntities[ codim ] = detail::makeSubEntities< Entity, codim >;
          } );
        cls.def( "subEntity", [ makeSubEntity ] ( const Entity &self, int i, int c ) {
            if( (c < codimension) || (c > dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( codimension ) + ", " + std::to_string( dimension ) + "])" );
            return makeSubEntity[ c ]( self, i );
          }, "index"_a, "codim"_a );
        cls.def( "subEntity", [ makeSubEntity ] ( const Entity &self, std::tuple< int, int > e ) {
            if( (std::get< 1 >( e ) < codimension) || (std::get< 1 >( e ) > dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( std::get< 1 >( e ) ) + " (must be in [" + std::to_string( codimension ) + ", " + std::to_string( dimension ) + "])" );
            return makeSubEntity[ std::get< 1 >( e ) ]( self, std::get< 0 >( e ) );
          } );
        cls.def( "subEntities", [ makeSubEntities ] ( const Entity &self, int c ) {
            if( (c < codimension) || (c > dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( codimension ) + ", " + std::to_string( dimension ) + "])" );
            return makeSubEntities[ c ]( self );
          }, "codim"_a );

        cls.def_property_readonly( "vertices", [] ( const Entity &self ) { return detail::makeSubEntities< Entity, dimension >( self ); } );
        if( mydimension < dimension )
        {
          cls.def_property_readonly( "edges", [] ( const Entity &self ) { return detail::makeSubEntities< Entity, (mydimension < dimension ? dimension-1 : dimension) >( self ); } );
          cls.def_property_readonly( "facets", [] ( const Entity &self ) { return detail::makeSubEntities< Entity, (mydimension < dimension ? mydimension+1 : dimension) >( self ); } );
        }
      }

      template< class Entity, class... options >
      inline static void registerExtendedEntityInterface ( pybind11::class_< Entity, options... > cls, PriorityTag< 0 > )
      {}

      template< class Entity, class... options >
      inline static void registerExtendedEntityInterface ( pybind11::class_< Entity, options... > cls )
      {
        return registerExtendedEntityInterface( cls, PriorityTag< 42 >() );
      }



      // registerGridEntity
      // ------------------

      template< class Entity, class... options >
      inline static void registerGridEntity ( pybind11::handle scope, pybind11::class_< Entity, options... > cls )
      {
        cls.def_property_readonly_static( "codimension", [] ( pybind11::object ) -> int { return Entity::codimension; } );
        cls.def_property_readonly_static( "dimension", [] ( pybind11::object ) -> int { return Entity::dimension; } );
        cls.def_property_readonly_static( "mydimension", [] ( pybind11::object ) -> int { return Entity::mydimension; } );

        cls.def_property_readonly( "geometry", [] ( const Entity &self ) { return self.geometry(); } );
        cls.def_property_readonly( "level", [] ( const Entity &self ) { return self.level(); } );
        cls.def_property_readonly( "type", [] ( const Entity &self ) { return self.type(); } );
        cls.def_property_readonly( "partitionType", [] ( const Entity &self ) { return self.partitionType(); } );
        cls.def_property_readonly( "referenceElement", [] ( const Entity &self ) { return referenceElement< double, Entity::dimension >( self.type() ); }, pybind11::keep_alive< 0, 1 >() );

        cls.def( pybind11::self == pybind11::self );
        cls.def( pybind11::self != pybind11::self );

        registerExtendedEntityInterface( cls );
      }

    } // namespace detail



    // registerGridEntity
    // ------------------

    template< class Grid, int codim >
    auto registerGridEntity ( pybind11::handle scope )
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;

      auto entry = insertClass<Entity>(scope, "Entity_" + std::to_string(codim),
         GenerateTypeName(MetaType<Grid>(), "template Codim<" + std::to_string(codim) + ">::Entity"));
      if (entry.second)
      {
        detail::registerGridEntity< Entity >( scope, entry.first );
        registerGridGeometry<Entity>( entry.first );
      }
      return entry.first;
    }

    // registerGridEntities
    // --------------------

    template< class Grid, int... codim >
    inline static auto registerGridEntities ( pybind11::handle scope, std::integer_sequence< int, codim... > )
    {
      return std::make_tuple( registerGridEntity< Grid, codim >( scope )... );
    }

    template< class Grid >
    auto registerGridEntities ( pybind11::handle scope )
    {
      return registerGridEntities< Grid >( scope, std::make_integer_sequence< int, Grid::dimension+1 >() );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_ENTITY_HH
