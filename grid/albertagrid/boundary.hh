// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_BOUNDARY_HH
#define DUNE_ALBERTA_BOUNDARY_HH

#include <dune/grid/albertagrid/albertaheader.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    template< int dim >
    struct BoundaryId;



#if DUNE_ALBERTA_VERSION >= 0x201
    template< int dim >
    struct BoundaryId
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_WALLS_MAX) );
        return elInfo->wall_bound[ face ];
      }
    };
#endif // #if DUNE_ALBERTA_VERSION >= 0x201



#if DUNE_ALBERTA_VERSION == 0x200
    template<>
    struct BoundaryId< 1 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_VERTICES_MAX) );
        return elInfo->vertex_bound[ face ];
      }
    };

    template<>
    struct BoundaryId< 2 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_EDGES_MAX) );
        return elInfo->edge_bound[ face ];
      }
    };

    template<>
    struct BoundaryId< 3 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_FACES_MAX) );
        return elInfo->face_bound[ face ];
      }
    };
#endif // #if DUNE_ALBERTA_VERSION == 0x200



#if DUNE_ALBERTA_VERSION < 0x200
#if DIM == 1
    template<>
    struct BoundaryId< 1 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_VERTICES) );
        return elInfo->bound[ face ];
      }
    };
#endif // #if DIM == 1

#if DIM == 2
    template<>
    struct BoundaryId< 2 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_EDGES) );
        return elInfo->boundary[ face ]->bound;
      }
    };
#endif // #if DIM == 2

#if DIM == 3
    template<>
    struct BoundaryId< 3 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int face )
      {
        assert( (face >= 0) && (face < N_FACES) );
        return elInfo->boundary[ face ]->bound;
      }
    };
#endif // #if DIM == 3
#endif // #if DUNE_ALBERTA_VERSION < 0x200



    template< int dim >
    inline int boundaryId ( ALBERTA EL_INFO *elInfo, int face )
    {
      return BoundaryId< dim >::get( elInfo, face );
    }

    inline bool isBoundary ( ALBERTA EL_INFO *elInfo, int face )
    {
      return (elInfo->neigh[ face ] == 0);
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_BOUNDARY_HH
