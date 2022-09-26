// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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

#include <dune/geometry/referenceelements.hh>

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
        if( !indexSet.contains( entity ) )
          pybind11::value_error( "Entity not contained in index set." );
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
          }, "entity"_a, "index"_a, "codim"_a,
          R"doc(
            get index of a subentity

            Numerical codes frequently require only the index of a subentity,
            not the subentity itself.
            This method directly provides the index.

            The subentity is described by a local index and a codimension.
            The exact interpretation of the local index is defined by the
            reference element.

            Calling this method is semantically equivalent to

            self.index(entity.subEntity(index, codim))

            Note: This method is implemented even if the grid does not implement
                  entities of the corresponding codimension.

            Args:
                entity:   entity containing the subentity
                index:    number of the subentity
                codim:    codimension of the subentity (wrt. the grid dimension)

            Returns:  Index of the subentity
          )doc" );
        cls.def( "subIndex", [] ( const IndexSet &self, const Entity &entity, std::tuple< int, int > e ) {
            return detail::indexSetSubIndex( self, entity, std::get< 0 >( e ), std::get< 1 >( e ) );
          }, "entity"_a, "subentity"_a,
          R"doc(
            Args:
                entity:       entity containing the subentity
                subentity:    (index, codim) tuple

            Returns:  Index of the subentity
          )doc" );

        cls.def( "subIndices", [] ( const IndexSet &self, const Entity &entity, int c ) {
            if( !self.contains( entity ) )
              pybind11::value_error( "Entity not contained in index set." );
            if( (c < Entity::codimension) || (c > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const int size = entity.subEntities( c );
            pybind11::tuple subIndices( size );
            for( int i = 0; i < size; ++i )
              subIndices[ i ] = pybind11::cast( self.subIndex( entity, i, c ) );
            return subIndices;
          }, "entity"_a, "codim"_a,
          R"doc(
            get indices of all subentities in a codimension

            Numerical codes frequently require all subindices of a codimension,
            e.g., all vertex indices.
            This convenience method provides them in a single call.

            Note: This method is implemented even if the grid does not implement
                  entities of the corresponding codimension.

            Args:
                entity:   entity containing the subentities
                codim:    codimension of the subentities (wrt. the grid dimension)

            Returns:  Tuple of indices, in the order given by the reference element
          )doc" );

      cls.def( "subIndices", [] ( const IndexSet &self, const Entity &entity, int i, int c, int cc ) {
            if( !self.contains( entity ) )
              pybind11::value_error( "Entity not contained in index set." );
            if( (c < Entity::codimension) || (c > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            if( (cc < c) || (cc > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( cc ) + " (must be in [" + std::to_string( c ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const auto reference = referenceElement< double, Entity::dimension >( entity.type() );
            const int size = reference.size( i, c, cc );
            pybind11::tuple subIndices( size );
            for( int ii = 0; ii < size; ++ii )
              subIndices[ ii ] = pybind11::cast( self.subIndex( entity, reference.subEntity(i, c, ii, cc), cc ) );
            return subIndices;
          }, "entity"_a, "index"_a, "codim"_a, "codim"_a,
          R"doc(
            get indices of all subentities in a codimension

            Numerical codes frequently require all subindices of a subentity of a codimension,
            e.g., all vertex indices of a edge.
            This convenience method provides them in a single call.

            Note: This method is implemented even if the grid does not implement
                  entities of the corresponding codimension.

            Args:
                entity:   entity containing the subentities
                index:    index of the subentity of entity
                codim:    codimension of the subentity of entity
                codim:    codimension of the subentities (wrt. the grid dimension)

            Returns:  Tuple of indices, in the order given by the reference element
          )doc" );

      cls.def( "subIndices", [] ( const IndexSet &self, const Entity &entity, std::pair< int, int > ic, int cc ) {
            int i = ic.first;
            int c = ic.second;
            if( !self.contains( entity ) )
              pybind11::value_error( "Entity not contained in index set." );
            if( (c < Entity::codimension) || (c > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            if( (cc < c) || (cc > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( cc ) + " (must be in [" + std::to_string( c ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const auto reference = referenceElement< double, Entity::dimension >( entity.type() );
            const int size = reference.size( i, c, cc );
            pybind11::tuple subIndices( size );
            for( int ii = 0; ii < size; ++ii )
              subIndices[ ii ] = pybind11::cast( self.subIndex( entity, reference.subEntity(i, c, ii, cc), cc ) );
            return subIndices;
          }, "entity"_a, "subentity"_a, "codim"_a,
          R"doc(
            get indices of all subentities in a codimension

            Numerical codes frequently require all subindices of a subentity of a codimension,
            e.g., all vertex indices of a edge.
            This convenience method provides them in a single call.

            Note: This method is implemented even if the grid does not implement
                  entities of the corresponding codimension.

            Args:
                entity:       entity containing the subentities
                subentity:    (index, codim) tuple
                codim:        codimension of the subentities (wrt. the grid dimension)

            Returns:  Tuple of indices, in the order given by the reference element
          )doc" );
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

    template< class IndexSet, class... options >
    inline static void registerGridViewIndexSet ( pybind11::handle scope, pybind11::class_< IndexSet, options... > cls )
    {
      using pybind11::operator""_a;

      pybind11::options opts;
      opts.disable_function_signatures();

      cls.doc() = R"doc(
          The index set maps all entities in a grid view to natural numbers
          (including 0).
          Denoting by E(g) the entities of geometry type g, this map has the
          following properties:
          - It is injective on each E(g).
          - The image of E(g) is consecutive and zero-starting, i.e.,
            { 0, ..., N-1 } for some integer N

          Index sets are used to assign user-defined data (e.g., degrees of
          freedom of a discretization) to entities in the grid view.
          For efficiency reasons, the preferred data structure for user data is
          an array.
          In order to access the data from associated to an entity, its index
          (wrt. an index set) is evaluated and used to access the  array element.

          Usually, the index set is not used directly. Instead, a mapper is used
          to compute the array index from the index provided by an index set.

          Note: The index assigned to an entity may change during grid
                modification (e.g., refinement of dynamic load balancing).
                The user is responsible for reorganizing the information stored
                in the external arrays.
                In order to do this, the concepts of id sets and persistent
                containers are provided.
        )doc";

      cls.def( "size", [] ( IndexSet &self, Dune::GeometryType type ) { return self.size( type ); } );
      cls.def( "size", [] ( IndexSet &self, int codim ) {
          if( (codim < 0) || (codim > IndexSet::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( IndexSet::dimension ) + "])" );
          return self.size( codim );
        } );

      cls.def( "types", [] ( const IndexSet &self, int codim ) {
          pybind11::list types;
          for( GeometryType type : self.types( codim ) )
            types.append( pybind11::cast( type ) );
          return types;
        }, "codim"_a,
        R"doc(
          obtain list of geometry types in the domain of the index set.

          Args:
              codim:    codimension to obtain geometry types for

          Returns:  List of all geometry types of given codimension in the
                    domain of the index set
        )doc" );

      Hybrid::forEach( std::make_integer_sequence< int, IndexSet::dimension+1 >(), [ &cls ] ( auto codim ) {
          typedef typename IndexSet::template Codim< codim >::Entity Entity;

          using pybind11::operator""_a;

          cls.def( "index", [] ( const IndexSet &self, const Entity &entity ) {
              if( !self.contains( entity ) )
                pybind11::value_error( "Entity not contained in index set." );
              return static_cast< int >( self.index( entity ) );
            }, "entity"_a,
            R"doc(
              get index of an entity

              Args:
                  entity:   entity to obtain index for

              Returns:  Index assigned to the entity
            )doc" );
          cls.def( "contains", [] ( const IndexSet &self, const Entity &entity ) {
              return self.contains( entity );
            }, "entity"_a,
            R"doc(
              check wheter an entity is in the domain of the index set

              Args:
                  entity:   entity to check for

              Returns: True, if the entity if in the domain, False otherwise
            )doc" );

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
      registerGridViewIndexSet( scope, cls );
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_INDEXSET_HH
