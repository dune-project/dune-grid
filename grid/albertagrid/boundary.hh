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

    template< int dim, int codim >
    struct BoundaryId;

    template< int dim >
    struct BoundaryId< dim, 0 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( subEntity == 0 );
        return 0;
      }
    };



#if DUNE_ALBERTA_VERSION >= 0x201
    struct VertexBoundaryId
    {
      static void clear ( ALBERTA EL_INFO *elInfo )
      {
        for( int i = 0; i < N_VERTICES_MAX; ++i )
          elInfo->vertex_bound[ i ] = 0;
      }

      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_VERTICES_MAX) );
        return elInfo->vertex_bound[ subEntity ];
      }
    };

#if DIM_MAX > 2
    struct EdgeBoundaryId
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_EDGES_MAX) );
        return elInfo->edge_bound[ subEntity ];
      }
    };
#endif // #if DIM_MAX > 2

    struct WallBoundaryId
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_WALLS_MAX) );
        return elInfo->wall_bound[ subEntity ];
      }
    };

    template<>
    struct BoundaryId< 1, 1 >
      : public VertexBoundaryId
    {};

    template<>
    struct BoundaryId< 2, 1 >
      : public WallBoundaryId
    {};

    template<>
    struct BoundaryId< 2, 2 >
      : public VertexBoundaryId
    {};

#if DIM_MAX > 2
    template<>
    struct BoundaryId< 3, 1 >
      : public WallBoundaryId
    {};

    template<>
    struct BoundaryId< 3, 2 >
      : public EdgeBoundaryId
    {};

    template<>
    struct BoundaryId< 3, 3 >
      : public VertexBoundaryId
    {};
#endif // #if DIM_MAX > 2
#endif // #if DUNE_ALBERTA_VERSION >= 0x201



#if DUNE_ALBERTA_VERSION == 0x200
    struct VertexBoundaryId
    {
      static void clear ( ALBERTA EL_INFO *elInfo )
      {
        for( int i = 0; i < N_VERTICES_MAX; ++i )
          elInfo->vertex_bound[ i ] = 0;
      }

      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_VERTICES_MAX) );
        return elInfo->vertex_bound[ subEntity ];
      }
    };

    struct EdgeBoundaryId
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_EDGES_MAX) );
        return elInfo->edge_bound[ subEntity ];
      }
    };

    struct FaceBoundaryId
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_FACES_MAX) );
        return elInfo->face_bound[ subEntity ];
      }
    };

    template<>
    struct BoundaryId< 1, 1 >
      : public VertexBoundaryId
    {};

    template<>
    struct BoundaryId< 2, 1 >
      : public EdgeBoundaryId
    {};

    template<>
    struct BoundaryId< 2, 2 >
      : public VertexBoundaryId
    {};

    template<>
    struct BoundaryId< 3, 1 >
      : public FaceBoundaryId
    {};

    template<>
    struct BoundaryId< 3, 2 >
      : public EdgeBoundaryId
    {};

    template<>
    struct BoundaryId< 3, 3 >
      : public VertexBoundaryId
    {};
#endif // #if DUNE_ALBERTA_VERSION == 0x200



#if DUNE_ALBERTA_VERSION < 0x200
    template<>
    struct BoundaryId< DIM, DIM >
    {
      static void clear ( ALBERTA EL_INFO *elInfo )
      {
        for( int i = 0; i < N_VERTICES; ++i )
          elInfo->bound[ i ] = 0;
      }

      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_VERTICES) );
        return elInfo->bound[ subEntity ];
      }
    };

#if DIM >= 2
    template<>
    struct BoundaryId< DIM, 1 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_FACES) );
        return elInfo->boundary[ subEntity ]->bound;
      }
    };
#endif // #if DIM >= 2

#if DIM >= 3
    template<>
    struct BoundaryId< DIM, 2 >
    {
      static int get ( ALBERTA EL_INFO *elInfo, int subEntity )
      {
        assert( (subEntity >= 0) && (subEntity < N_EDGES) );
        return elInfo->boundary[ N_FACES + subEntity ]->bound;
      }
    };
#endif // #if DIM >= 3
#endif // #if DUNE_ALBERTA_VERSION < 0x200



    template< int dim, int codim >
    inline int boundaryId ( ALBERTA EL_INFO *elInfo, int subEntity )
    {
      return BoundaryId< dim, codim >::get( elInfo, subEntity );
    }

    inline bool isBoundary ( ALBERTA EL_INFO *elInfo, int face )
    {
      return (elInfo->neigh[ face ] == 0);
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_BOUNDARY_HH
