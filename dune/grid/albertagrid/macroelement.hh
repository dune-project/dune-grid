// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MACROELEMENT_HH
#define DUNE_ALBERTA_MACROELEMENT_HH

#include <dune/grid/albertagrid/misc.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // MacroElement
    // ------------

    template< int dim >
    struct MacroElement
      : public ALBERTA MACRO_EL
    {
      const GlobalVector &coordinate ( const int vertex ) const;

      int boundaryId ( const int face ) const;
      bool isBoundary ( const int face ) const;
      const MacroElement< dim > *neighbor ( const int face ) const;
    };


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline const GlobalVector &MacroElement< dim >::coordinate ( const int vertex ) const
    {
      assert( (vertex >= 0) && (vertex < N_VERTICES_MAX) );
      return *coord[ vertex ];
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300


#if DUNE_ALBERTA_VERSION < 0x300
    template< int dim >
    inline const GlobalVector &MacroElement< dim >::coordinate ( const int vertex ) const
    {
      assert( (vertex >= 0) && (vertex < N_VERTICES_MAX) );
      return *((const GlobalVector *)coord[ vertex ]);
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x300


    template< int dim >
    inline bool MacroElement< dim >::isBoundary ( const int face ) const
    {
      return (boundaryId( face ) != InteriorBoundary);
    }


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline int MacroElement< dim >::boundaryId ( const int face ) const
    {
      return wall_bound[ face ];
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION < 0x300
    template< int dim >
    inline int MacroElement< dim >::boundaryId ( const int face ) const
    {
      switch( dim )
      {
      case 1 :
        assert( (face >= 0) && (face < N_VERTICES_MAX) );
        return vertex_bound[ face ];
      case 2 :
        assert( (face >= 0) && (face < N_EDGES_MAX) );
        return edge_bound[ face ];
      case 3 :
        assert( (face >= 0) && (face < N_FACES_MAX) );
        return face_bound[ face ];
      }
    }
#endif // #if DUNE_ALBERTA_VERSION < 0x300


    template< int dim >
    const MacroElement< dim > *MacroElement< dim >::neighbor ( const int face ) const
    {
      assert( (face >= 0) && (face < N_NEIGH_MAX) );
      return static_cast< const MacroElement * >( neigh[ face ] );
    }

  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_MACROELEMENT_HH
