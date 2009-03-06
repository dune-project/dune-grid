// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_REFERENCETOPO_HH
#define DUNE_ALBERTAGRID_REFERENCETOPO_HH

#if HAVE_ALBERTA

#include <cassert>

#include <dune/grid/albertagrid/albertaheader.hh>

#ifdef __ALBERTApp__
namespace Albert {
#endif

namespace AlbertHelp {

  // NOTE: Vertex numbering in ALBERTA is the same as in Dune
  // therefore no map is provided for that

  // faces in 2d (i.e. triangle edges )
  // which vertices belong to which face
  static const int localTriangleFaceNumber [3][2] = { {1,2} , {2,0} , {0,1} };

  // edges in 3d
  // local numbers of vertices belonging to one edge
  // according to Alberta reference element which is for edges different to
  // the Dune reference simplex in 3d ,see Alberta doc page 105
  static const int localEdgeNumber [6][2] =
  {
    {0,1} , {0,2} , {0,3} , // first three vertices like in 2d for faces(edges)
    {1,2} , {1,3} , {2,3} // then all with the last vertex
  };

  // see Albert Doc page 12 for reference element
  // if we look from outside, then face numbering must be clockwise
  // see below that the vertex numbers for each face are the same
  // but in Dune reference element not clockwise ,but this we need for
  // calculation the outer normal, see calcOuterNormal in albertagrid.cc
  static const int tetraFace_0[3] = {1,3,2};
  static const int tetraFace_1[3] = {0,2,3};
  static const int tetraFace_2[3] = {0,3,1};
  static const int tetraFace_3[3] = {0,1,2};

  // local vertex numbering of faces in the Dune refrence simplex in 3d
  static const int localDuneTetraFaceNumber[4][3] =
  { {1,2,3}, // face 0
    {0,3,2}, // face 1
    {0,1,3}, // face 2
    {0,2,1} // face 3
  };

  static const int *const localAlbertaFaceNumber[ 4 ]
    = { tetraFace_0, tetraFace_1, tetraFace_2 , tetraFace_3 };

  //****************************************************************
  //
  //  specialization of mapVertices
  //  see referencetopo.hh
  //
  //****************************************************************
  template <int md, int cd>
  struct MapVertices
  {
    static int mapVertices ( int subEntity, int i )
    {
      assert( subEntity == 0 );
      return i;
    }
  };

  // faces in 2d
  template <>
  struct MapVertices<1,2>
  {
    static int mapVertices ( int subEntity, int i )
    {
      assert( (subEntity >= 0) && (subEntity < 3) );
      assert( (i >= 0) && (i < 2) );
      return ALBERTA AlbertHelp :: localTriangleFaceNumber[ subEntity ][ i ];
    }
  };

  // faces in 3d
  template <>
  struct MapVertices<2,3>
  {
    static int mapVertices ( int subEntity, int i )
    {
      assert( (subEntity >= 0) && (subEntity < 4) );
      assert( (i >= 0) && (i < 3) );
      return ALBERTA AlbertHelp :: localDuneTetraFaceNumber[ subEntity ][ i ];
    }
  };

  // edges in 3d
  template <>
  struct MapVertices<1,3>
  {
    static int mapVertices ( int subEntity, int i )
    {
      assert( (subEntity >= 0) && (subEntity < 6) );
      assert( (i >= 0) && (i < 2) );
      return ALBERTA AlbertHelp :: localEdgeNumber[ subEntity ][ i ];
    }
  };

  // vertices in 2d and 3d
  template <int cd>
  struct MapVertices<0,cd>
  {
    static int mapVertices ( int subEntity, int i )
    {
      assert( i == 0 );
      return subEntity;
    }
  };


} // end namespace AlbertHelp

#ifdef __ALBERTApp__
} // end namespace Albert
#endif

#endif // HAVE_ALBERTA

#endif
