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
    for( const auto vertex : vertices( grid.leafGridView() ) )
    {
      std::size_t idx = factory.insertionIndex( vertex );
      Vertex v = mesh.vertices[ idx ];
      if( (v - vertex.geometry().center() ).two_norm() > 1e-8 )
        DUNE_THROW( GridError, "GridFactory error, Vertex insertion Index wrong!" );
    }

    // check element insertion index
    std::vector< unsigned int > insertIndices, indices;
    for( const auto element : elements( grid.leafGridView() ) )
    {
      std::size_t idx = factory.insertionIndex( element );
      unsigned int numSubeEntitites = element.subEntities( Grid::dimension );

      insertIndices = mesh.elements[ idx ].second;
      std::sort( insertIndices.begin(), insertIndices.end() );

      if( numSubeEntitites != insertIndices.size() )
        DUNE_THROW( GridError, "GridFactory error, wrong number of subEntities inserted!" );

      indices.clear();
      for( unsigned int i = 0; i < numSubeEntitites; ++i )
        indices.push_back( factory.insertionIndex( element.template subEntity< Grid::dimension >( i ) ) );
      std::sort( indices.begin(), indices.end() );

      if( !std::equal( indices.begin(), indices.end(), insertIndices.begin() ) )
        DUNE_THROW( GridError, "GridFactory error, Element insertion index wrong!" );
    }
  }

  template< class Grid, class Mesh >
  void checkGridFactory ( const Mesh &mesh )
  {
    checkGridFactory< Grid >( mesh, [] ( const typename Mesh::Vertex &v ) { return v; } );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKGRIDFACTORY_HH
