// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/alugrid/2d/geoinfather.hh>

// coordinates for the children of split4 quadrilaterals
static const Dune::alu2d_ctype quadrilateral[ 4 ][ 4 ][ 2 ]
  = { { { 0.0, 0.0 }, { 0.5, 0.0 }, { 0.0, 0.5 }, { 0.5, 0.5 } },
      { { 0.5, 0.0 }, { 1.0, 0.0 }, { 0.5, 0.5 }, { 1.0, 0.5 } },
      { { 0.5, 0.5 }, { 1.0, 0.5 }, { 0.5, 1.0 }, { 1.0, 1.0 } },
      { { 0.0, 0.5 }, { 0.5, 0.5 }, { 0.0, 1.0 }, { 0.5, 1.0 } } };

// coordinates for the children of split4 triangles
static const Dune::alu2d_ctype triangle_nonconforming[ 4 ][ 3 ][ 2 ]
  = { { { 0.0, 0.0 }, { 0.5, 0.0 }, { 0.0, 0.5 } },
      { { 0.5, 0.0 }, { 1.0, 0.0 }, { 0.5, 0.5 } },
      { { 0.0, 0.5 }, { 0.5, 0.5 }, { 0.0, 1.0 } },
      { { 0.5, 0.5 }, { 0.0, 0.5 }, { 0.5, 0.0 } } };

// coordinates for the children of split2 triangles
static const Dune::alu2d_ctype triangle_conforming[ 2 ][ 3 ][ 2 ]
  = { { { 0.5, 0.5 }, { 0.0, 0.0 }, { 1.0, 0.0 } },
      { { 0.5, 0.5 }, { 0.0, 1.0 }, { 0.0, 0.0 } } };

namespace Dune
{

  void getALU2DGeometryInFather ( bool conforming, int child, array< FieldVector< alu2d_ctype, 2 >, 3 > &corners )
  {
    if( conforming )
    {
      assert( (child >= 0) && (child < 2) );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 2; ++j )
          corners[ i ][ j ] = triangle_conforming[ child ][ i ][ j ];
      }
    }
    else
    {
      assert( (child >= 0) && (child < 4) );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 2; ++j )
          corners[ i ][ j ] = triangle_nonconforming[ child ][ i ][ j ];
      }
    }
  }

  void getALU2DGeometryInFather ( int child, array< FieldVector< alu2d_ctype, 2 >, 4 > &corners )
  {
    assert( (child >= 0) && (child < 4) );
    for( int i = 0; i < 4; ++i )
    {
      for( int j = 0; j < 2; ++j )
        corners[ i ][ j ] = quadrilateral[ child ][ i ][ j ];
    }
  }

} // namespace Dune
