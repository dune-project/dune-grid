// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_TEST_CHECKGRIDFACTORY_HH
#define DUNE_GRID_TEST_CHECKGRIDFACTORY_HH

#include <algorithm>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>

namespace Dune
{

  // checkGridFactory
  // ----------------

  template< class Grid, class Mesh, class Projection >
  void checkGridFactory ( const Mesh &mesh, Projection &&projection )
  {
    GridFactory< Grid > factory;

    // create grid from mesh
    typedef FieldVector< typename Grid::ctype, Grid::dimensionworld > Vertex;

    mesh.addToGridFactory( factory, projection );

    std::unique_ptr< Grid > gridptr( factory.createGrid() );
    Grid &grid = *gridptr;

    // check insertion indices

    // check vertex insertion index
    for( const auto &vertex : vertices( grid.leafGridView() ) )
    {
      std::size_t idx = factory.insertionIndex( vertex );
      Vertex v = projection( mesh.vertices[ idx ] );
      if( (v - vertex.geometry().center() ).two_norm() > 1e-8 )
        DUNE_THROW( GridError, "GridFactory error, Vertex insertion Index wrong!" );
    }

    // check element insertion index
    std::vector< unsigned int > indices;
    for( const auto &element : elements( grid.leafGridView() ) )
    {
      std::size_t idx = factory.insertionIndex( element );
      unsigned int numSubEntitites = element.subEntities( Grid::dimension );

      if( numSubEntitites != mesh.elements[ idx ].second.size() )
        DUNE_THROW( GridError, "GridFactory error, wrong number of subEntities inserted!" );

      indices.clear();
      for( unsigned int i = 0; i < numSubEntitites; ++i )
        indices.push_back( factory.insertionIndex( element.template subEntity< Grid::dimension >( i ) ) );

      if( !std::is_permutation( indices.begin(), indices.end(), mesh.elements[ idx ].second.begin() ) )
        DUNE_THROW( GridError, "GridFactory error, Element insertion index wrong!" );
    }

    // check boundary segment index
    typedef std::pair< unsigned int, std::vector< unsigned int > > BoundarySegementPair;
    std::vector< BoundarySegementPair > bndInsIndex;

    for( const auto &entity : elements( grid.leafGridView() ) )
    {
      auto refelement = referenceElement< typename Grid::ctype, Grid::dimension >( entity.type() );

      for( const auto &intersection : intersections( grid.leafGridView(), entity ) )
        if( factory.wasInserted( intersection ) )
        {
          std::vector< unsigned int > vertices;
          const int numSubEntitites = refelement.size( intersection.indexInInside(), 1, Grid::dimension );
          for( int i = 0; i < numSubEntitites; ++i )
            vertices.push_back(
              factory.insertionIndex(
                entity.template subEntity< Grid::dimension >( refelement.subEntity( intersection.indexInInside(), 1, i, Grid::dimension ) )
                ) );
          bndInsIndex.emplace_back( factory.insertionIndex( intersection ), vertices );
        }
    }

    auto compare = [] ( const BoundarySegementPair &v, const BoundarySegementPair &w ) { return v.first == w.first; };

    bndInsIndex.resize( std::distance( bndInsIndex.begin(), std::unique( bndInsIndex.begin(), bndInsIndex.end(), compare ) ) );

    if( bndInsIndex.size() != mesh.boundaries.size() )
      DUNE_THROW( GridError, "GridFactory error, wrong number of boundary segments inserted." );

    for( const BoundarySegementPair &pair : bndInsIndex )
      if( !std::is_permutation( pair.second.begin(), pair.second.end(), mesh.boundaries[ pair.first ].begin() ) )
        DUNE_THROW( GridError, "GridFactory error, insertion index for boundary segment wrong." );
  }


  template< class Grid, class Mesh >
  void checkGridFactory ( const Mesh &mesh )
  {
    checkGridFactory< Grid >( mesh, [] ( const typename Mesh::Vertex &v ) { return v; } );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKGRIDFACTORY_HH
