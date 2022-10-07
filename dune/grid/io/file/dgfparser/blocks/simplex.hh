// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_SIMPLEXBLOCK_HH
#define DUNE_DGF_SIMPLEXBLOCK_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {
    // SimplexBlock
    // ------------

    class SimplexBlock
      : public BasicBlock
    {
      unsigned int nofvtx;
      int vtxoffset;
      int dimgrid;
      bool goodline;                 // active line describes a vertex
      int nofparams;                 // nof parameters

    public:
      SimplexBlock ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid );

      int get ( std :: vector< std :: vector< unsigned int > > &simplex,
                std :: vector< std :: vector< double > > &params,
                int &nofp );

      // cubes -> simplex
      static int
      cube2simplex ( std :: vector< std :: vector< double > > &vtx,
                     std :: vector< std :: vector< unsigned int > > &elements,
                     std :: vector< std :: vector< double > > &params );

      // some information
      bool ok ()
      {
        return goodline;
      }

      int nofsimplex ()
      {
        return noflines();
      }

    private:
      // get the dimension of the grid
      int getDimGrid ();
      // get next simplex
      bool next ( std :: vector< unsigned int > &simplex,
                  std :: vector< double > &param );
    };

  } // end namespace dgf

} // end namespace Dune

#endif
