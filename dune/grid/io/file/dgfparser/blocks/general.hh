// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_DGFPARSER_BLOCKS_GENERAL_HH
#define DUNE_GRID_IO_FILE_DGFPARSER_BLOCKS_GENERAL_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {

    // GeneralBlock
    // ---------

    class GeneralBlock
      : public BasicBlock
    {
      unsigned int nofvtx;
      int dimgrid;
      bool goodline;        // active line describes a vertex
      std :: vector< unsigned int > map; // active vertex
      int nofparams;
      int vtxoffset;

    public:
      GeneralBlock ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid );

      int get ( std :: vector< std :: vector< unsigned int> > &simplex,
                std :: vector< std :: vector< double > > &params,
                int &nofp );

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

#endif // #ifndef DUNE_GRID_IO_FILE_DGFPARSER_BLOCKS_GENERAL_HH
