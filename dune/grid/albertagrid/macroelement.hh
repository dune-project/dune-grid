// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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


    template< int dim >
    inline const GlobalVector &MacroElement< dim >::coordinate ( const int vertex ) const
    {
      assert( (vertex >= 0) && (vertex < N_VERTICES_MAX) );
      return *coord[ vertex ];
    }


    template< int dim >
    inline bool MacroElement< dim >::isBoundary ( const int face ) const
    {
      return (boundaryId( face ) != InteriorBoundary);
    }


    template< int dim >
    inline int MacroElement< dim >::boundaryId ( const int face ) const
    {
      return wall_bound[ face ];
    }


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
