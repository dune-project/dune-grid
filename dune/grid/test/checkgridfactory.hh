#ifndef DUNE_GRID_TEST_CHECKGRIDFACTORY_HH
#define DUNE_GRID_TEST_CHECKGRIDFACTORY_HH

#include <algorithm>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune
{

  // checkGridFactory
  // ----------------

  template< class Grid, class Mesh >
  void checkGridFactory ( const Mesh &mesh )
  {
    GridFactory< Grid > factory;

    // create grid from mesh
    typedef FieldVector< typename Grid::ctype, Grid::dimensionworld > Vertex;

    for( Vertex vertex : mesh.vertices() )
      factory.insertVertex( vertex );

    for( const auto &element : mesh.elements() )
      factory.insertElement( element.type(), element.vertices() );

    std::unique_ptr< Grid > gridptr( factory.createGrid() );
    Grid &grid = *gridptr;

    // check insertion indices

    // check vertex insertion index
    for( const auto vertex : vertices( grid.leafGridView() ) )
    {
      std::size_t idx = factory.insertionIndex( vertex );
      Vertex v = mesh.vertices()[ idx ];
      if( (v - vertex.geometry().center() ).two_norm() > 1e-8 )
        DUNE_THROW( GridError, "GridFactor error, Vertex insertion Index wrong!" );
    }

    // check element insertion index
    std::vector< unsigned int > insertIndices, indices;
    for( const auto element : elements( grid.leafGridView() ) )
    {
      indices.clear();
      insertIndices.clear();

      std::size_t idx = factory.insertionIndex( element );
      unsigned int numSubeEntitites = element.subEntities( Grid::dimension );

      if( numSubeEntitites != mesh.elements()[ idx ].vertices().size() )
        DUNE_THROW( GridError, "GridFactor error, wrong number of subEntities inserted!" );

      for( unsigned int i = 0; i < numSubeEntitites; ++i )
      {
        indices.push_back( factory.insertionIndex( element.template subEntity< Grid::dimension >( i ) ) );
        insertIndices.push_back( mesh.elements()[ idx ].vertices()[ i ] );
      }

      std::sort( insertIndices.begin(), insertIndices.end() );
      std::sort( indices.begin(), indices.end() );

      for( unsigned int i = 0; i < indices.size(); ++i )
        if( indices[ i ] != insertIndices[ i ] )
          DUNE_THROW( GridError, "GridFactor error, Element insertion index wrong!" );
    }

  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKGRIDFACTORY_HH
