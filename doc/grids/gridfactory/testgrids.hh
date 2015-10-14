#ifndef DOC_GRIDS_GRIDFACTORY_TESTGRIDS_HH
#define DOC_GRIDS_GRIDFACTORY_TESTGRIDS_HH

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridfactory.hh>

namespace Dune
{

  // TestGrid
  // --------

  template< int dimworld >
  struct TestGrid
  {
    typedef FieldVector< double, dimworld > Vertex;
    typedef std::pair< GeometryType, std::vector< unsigned int > > Element;
    typedef std::vector< unsigned int > BoundarySegment;

    TestGrid ( std::vector< Vertex > v, std::vector< Element > e, std::vector< BoundarySegment > b )
      : vertices( std::move( v ) ), elements( std::move( e ) ), boundaries( std::move( b ) )
    {}

    template< class Grid, class Projection >
    void addToGridFactory ( Dune::GridFactory< Grid > &factory, Projection &&projection ) const
    {
      for( const Vertex &v : vertices )
        factory.insertVertex( projection( v ) );
      for( const Element &e : elements )
        factory.insertElement( e.first, e.second );
      for( const BoundarySegment &b : boundaries )
        factory.insertBoundarySegment( b );
    }

    template< class Grid >
    void addToGridFactory ( Dune::GridFactory< Grid > &factory ) const
    {
      addToGridFactory( factory, [] ( const Vertex &v ) { return v; } );
    }

    std::vector< Vertex > vertices;
    std::vector< Element > elements;
    std::vector< BoundarySegment > boundaries;
  };


  namespace TestGrids
  {

    static const GeometryType triangle( GeometryType::simplex, 2 );
    static const GeometryType quadrilateral( GeometryType::cube, 2 );

    static const GeometryType tetrahedron( GeometryType::simplex, 3 );
    static const GeometryType hexahedron( GeometryType::cube, 3 );

    // line grids

    static const TestGrid< 1 > unitLine = {
        { { 0.0 }, { 1.0 } },
        { { line, { 0, 1 } } },
        { { 0 }, { 1 } }
      };


    // cube grids
    static const TestGrid< 2 > unitSquare = {
        { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } },
        { { quadrilateral, { 0, 1, 2, 3 } } },
        { { 0, 1 }, { 1, 3 }, { 3, 2 }, { 2, 0 } }
      };

    static const TestGrid< 2 > unitCube = {
        {
          { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
        },
        { { hexahedron, { 0, 1, 2, 3, 4, 5, 6, 7 } } },
        {
          { 0, 1, 2, 3 }, { 1, 3, 7, 5 },
          { 0, 4, 5, 1 }, { 0, 4, 6, 2 },
          { 2, 3, 7, 6 }, { 4, 5, 7, 6 }
        }
      };

    static const TestGrid< 2 > kuhn2d = {
        { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } },
        { { triangle, { 0, 1, 3 } }, { triangle, { 0, 2, 3 } } },
        { { 0, 1 }, { 1, 3 }, { 3, 2 }, { 2, 0 } }
      };

    static const TestGrid< 3 > kuhn3d = {
        {
          { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
        },
        {
          { tetrahedron, { 0, 1, 3, 7 } }, { tetrahedron, { 0, 1, 5, 7 } },
          { tetrahedron, { 0, 2, 3, 7 } }, { tetrahedron, { 0, 2, 6, 7 } },
          { tetrahedron, { 0, 4, 5, 7 } }, { tetrahedron, { 0, 4, 6, 7 } }
        },
        {
          { 0, 1, 3 }, { 1, 3, 7 },
          { 0, 1, 5 }, { 1, 5, 7 },
          { 0, 2, 3 }, { 2, 3, 7 },
          { 0, 2, 6 }, { 2, 6, 7 },
          { 0, 4, 5 }, { 4, 5, 7 },
          { 0, 4, 6 }, { 4, 6, 7 }
        }
      };

  } // namespace TestGrids

} // namespace Dune

#endif // #ifndef DOC_GRIDS_GRIDFACTORY_TESTGRIDS_HH
