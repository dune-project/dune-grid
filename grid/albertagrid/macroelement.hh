// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MACROELEMENT_HH
#define DUNE_ALBERTA_MACROELEMENT_HH

#include <dune/grid/albertagrid/misc.hh>

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
      bool isBoundary ( const int face ) const;
      const MacroElement< dim > *neighbor ( const int face ) const;
    };


#if DUNE_ALBERTA_VERSION >= 0x300
    template< int dim >
    inline bool MacroElement< dim >::isBoundary ( const int face ) const
    {
      const int id = wall_bound[ face ];
      return (id != 0);
    }
#endif // #if DUNE_ALBERTA_VERSION >= 0x300

#if DUNE_ALBERTA_VERSION == 0x200
    template<>
    inline bool MacroElement< 1 >::isBoundary ( const int face ) const
    {
      assert( (face >= 0) && (face < N_VERTICES_MAX) );
      return vertex_bound[ face ];
    }

    template<>
    inline bool MacroElement< 2 >::isBoundary ( const int face ) const
    {
      assert( (face >= 0) && (face < N_EDGES_MAX) );
      return edge_bound[ face ];
    }

    template<>
    inline bool MacroElement< 3 >::isBoundary ( const int face ) const
    {
      assert( (face >= 0) && (face < N_FACES_MAX) );
      return face_bound[ face ];
    }
#endif // #if DUNE_ALBERTA_VERSION == 0x200


    template< int dim >
    const MacroElement< dim > *MacroElement< dim >::neighbor ( const int face ) const
    {
      assert( (face >= 0) && (face < N_FACES_MAX) );
      return static_cast< const MacroElement * >( neigh[ face ] );
    }


  }

}


#endif // #ifndef DUNE_ALBERTA_MACROELEMENT_HH
