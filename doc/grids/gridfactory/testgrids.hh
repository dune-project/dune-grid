// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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

    using namespace GeometryTypes;

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

    static const TestGrid< 3 > unitCube = {
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

    // simplex grids
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

    // hybrid grids
    static const TestGrid< 2 > hybrid2d = {
        {
          {0, 0}, {0.5, 0}, {0.5, 0.5}, {0, 0.5},
          {0.25, 0}, {0.5, 0.25}, {0.25, 0.5}, {0, 0.25},
          {0.25, 0.25}, {1, 0}, {1, 0.5}, {0.75, 0.25},
          {1, 1}, {0.5, 1}, {0, 1}, {0.25, 0.75}
        },
        {
          { triangle, { 9, 10, 11} }, { triangle, { 15, 13, 14 } },
          { quadrilateral, {0, 4, 7, 8} }, { quadrilateral, {4, 1, 8, 5} },
          { quadrilateral, {8, 5, 6, 2} }, { quadrilateral,{7, 8, 3, 6} },
          { quadrilateral, {1, 9, 5, 11} }, { quadrilateral, {5, 11, 2, 10} },
          { quadrilateral, {2, 10, 13, 12} }, { quadrilateral, {3, 6, 14, 15} },
          { quadrilateral, {6, 2, 15, 13} }
        },
        {}
      };

   static const TestGrid< 3 > hybrid3d = {
        {
          {0, 0, 0}, {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0},
          {0, 0, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0.5}, {0, 0.5, 0.5},
          {0.25, 0, 0}, {0.5, 0.25, 0}, {0.25, 0.5, 0}, {0, 0.25, 0},
          {0, 0, 0.25}, {0.5, 0, 0.25}, {0.5, 0.5, 0.25}, {0, 0.5, 0.25},
          {0.25, 0, 0.5}, {0.5, 0.25, 0.5}, {0.25, 0.5, 0.5}, {0, 0.25, 0.5},
          {0.25, 0.25, 0}, {0.25, 0, 0.25}, {0.5, 0.25, 0.25}, {0.25, 0.5, 0.25},
          {0, 0.25, 0.25}, {0.25, 0.25, 0.5}, {0.25, 0.25, 0.25}, {1, 0, 0},
          {1, 0.5, 0}, {1, 0, 0.5}, {1, 0.5, 0.5}, {0.75, 0.25, 0.25},
          {1, 1, 0}, {0.5, 1, 0}, {1, 1, 0.5}, {0.5, 1, 0.5},
          {0.75, 0.75, 0.25}, {0, 1, 0}, {0, 1, 0.5}, {0.25, 0.75, 0.25},
          {0, 0, 1}, {0.5, 0, 1}, {0.5, 0.5, 1}, {0, 0.5, 1},
          {0.25, 0.25, 0.75}, {1, 0, 1}, {1, 0.5, 1}, {0.75, 0.25, 0.75},
          {1, 1, 1}, {0.5, 1, 1}, {0, 1, 1}, {0.25, 0.75, 0.75},
          {1.5, 0, 0}, {1.5, 0.5, 0}, {1.5, 1, 0}, {1.5, 0, 0.5},
          {1.5, 0.5, 0.5}, {1.5, 1, 0.5}, {1.5, 0, 1}, {1.5, 0.5, 1},
          {1.5, 1, 1}
        },
        {
          { tetrahedron,{10, 29, 3, 32} }, { tetrahedron,{10, 2, 28, 32}},
          { tetrahedron,{10, 28, 29, 32}}, { tetrahedron,{14, 28, 2, 32}},
          { tetrahedron,{14, 6, 30, 32}}, { tetrahedron,{14, 30, 28, 32}},
          { tetrahedron,{15, 31, 7, 32}}, { tetrahedron,{15, 3, 29, 32}},
          { tetrahedron,{15, 29, 31, 32}}, { tetrahedron,{18, 30, 6, 32}},
          { tetrahedron,{18, 7, 31, 32}}, { tetrahedron,{18, 31, 30, 32}},
          { tetrahedron,{15, 29, 3, 37}}, { tetrahedron,{15, 7, 31, 37}},
          { tetrahedron,{15, 31, 29, 37}}, { tetrahedron,{15, 36, 7, 37}},
          { tetrahedron,{15, 3, 34, 37}}, { tetrahedron,{15, 34, 36, 37}},
          { tetrahedron,{11, 38, 4, 40}}, { tetrahedron,{11, 3, 34, 40}},
          { tetrahedron,{11, 34, 38, 40}}, { tetrahedron,{15, 34, 3, 40}},
          { tetrahedron,{15, 7, 36, 40}}, { tetrahedron,{15, 36, 34, 40}},
          { tetrahedron,{16, 39, 8, 40}}, { tetrahedron,{16, 4, 38, 40}},
          { tetrahedron,{16, 38, 39, 40}}, { tetrahedron,{19, 36, 7, 40}},
          { tetrahedron,{19, 8, 39, 40}}, { tetrahedron,{19, 39, 36, 40}},
          { tetrahedron,{17, 42, 6, 45}}, { tetrahedron,{17, 5, 41, 45}},
          { tetrahedron,{17, 41, 42, 45}}, { tetrahedron,{18, 43, 7, 45}},
          { tetrahedron,{18, 6, 42, 45}}, { tetrahedron,{18, 42, 43, 45}},
          { tetrahedron,{19, 44, 8, 45}}, { tetrahedron,{19, 7, 43, 45}},
          { tetrahedron,{19, 43, 44, 45}}, { tetrahedron,{20, 41, 5, 45}},
          { tetrahedron,{20, 8, 44, 45}}, { tetrahedron,{20, 44, 41, 45}},
          { tetrahedron,{18, 31, 7, 48}}, { tetrahedron,{18, 6, 30, 48}},
          { tetrahedron,{18, 30, 31, 48}}, { tetrahedron,{18, 42, 6, 48}},
          { tetrahedron,{18, 7, 43, 48}}, { tetrahedron,{18, 43, 42, 48}},
          { tetrahedron,{19, 39, 8, 52}}, { tetrahedron,{19, 7, 36, 52}},
          { tetrahedron,{19, 36, 39, 52}}, { tetrahedron,{19, 43, 7, 52}},
          { tetrahedron,{19, 8, 44, 52}}, { tetrahedron,{19, 44, 43, 52}},
          { pyramid, {28, 30, 29, 31, 32}}, { pyramid, {10, 23, 2,  14, 32}},
          { pyramid, {14, 23, 6,  18, 32}}, { pyramid, {18, 23, 7,  15, 32}},
          { pyramid, {15, 23, 3,  10, 32}}, { pyramid, {3,  29, 34, 33, 37}},
          { pyramid, {29, 31, 33, 35, 37}}, { pyramid, {33, 35, 34, 36, 37}},
          { pyramid, {7,  36, 31, 35, 37}}, { pyramid, {11, 24, 3,  15, 40}},
          { pyramid, {15, 24, 7,  19, 40}}, { pyramid, {19, 24, 8,  16, 40}},
          { pyramid, {16, 24, 4,  11, 40}}, { pyramid, {34, 36, 38, 39, 40}},
          { pyramid, {20, 26, 8,  19, 45}}, { pyramid, {19, 26, 7,  18, 45}},
          { pyramid, {18, 26, 6,  17, 45}}, { pyramid, {17, 26, 5,  20, 45}},
          { pyramid, {41, 44, 42, 43, 45}}, { pyramid, {6,  42, 30, 46, 48}},
          { pyramid, {30, 46, 31, 47, 48}}, { pyramid, {31, 47, 7,  43, 48}},
          { pyramid, {42, 43, 46, 47, 48}}, { pyramid, {7,  43, 36, 50, 52}},
          { pyramid, {36, 50, 39, 51, 52}}, { pyramid, {39, 51, 8,  44, 52}},
          { pyramid, {44, 51, 43, 50, 52}},
          { prism, {28, 53, 29, 30, 56, 31}}, { prism, {53, 54, 29, 56, 57, 31}},
          { prism, {30, 56, 31, 46, 59, 47}}, { prism, {56, 57, 31, 59, 60, 47}},
          { prism, {29, 54, 33, 31, 57, 35}}, { prism, {54, 55, 33, 57, 58, 35}},
          { prism, {31, 57, 35, 47, 60, 49}}, { prism, {57, 58, 35, 60, 61, 49}},
          { hexahedron, {1, 9, 12, 21, 13, 22, 25, 27}}, { hexahedron, {9, 2, 21, 10, 22, 14, 27, 23}},
          { hexahedron, {21, 10, 11, 3, 27, 23, 24, 15}}, { hexahedron, {12, 21, 4, 11, 25, 27, 16, 24}},
          { hexahedron, {13, 22, 25, 27, 5, 17, 20, 26}}, { hexahedron, {22, 14, 27, 23, 17, 6, 26, 18}},
          { hexahedron, {27, 23, 24, 15, 26, 18, 19, 7}}, { hexahedron, {25, 27, 16, 24, 20, 26, 8, 19}},
          { hexahedron, {7, 31, 36, 35, 43, 47, 50, 49}}
        },
        {}
      };

  } // namespace TestGrids

} // namespace Dune

#endif // #ifndef DOC_GRIDS_GRIDFACTORY_TESTGRIDS_HH
