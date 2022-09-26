// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_VERTEXBLOCK_HH
#define DUNE_DGF_VERTEXBLOCK_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {

    class VertexBlock
      : public BasicBlock
    {
      int dimvertex;         // the dimension of the vertices (determined from DGF file)
      int dimworld;          // the dimension of the world (either dimvertex or given by user)
      bool goodline;         // active line describes a vertex
      int vtxoffset;
      int nofParam;

    public:
      // initialize vertex block
      VertexBlock ( std :: istream &in, int &pdimworld );

      int get ( std :: vector< std :: vector< double > > &vtx,
                std :: vector< std :: vector< double > > &param,
                int &nofp );

      // some information
      bool ok () const
      {
        return goodline;
      }

      int offset () const
      {
        return vtxoffset;
      }

    private:
      // get dimworld
      int getDimWorld ();

      // get next vertex
      bool next ( std :: vector< double > &point, std :: vector< double > &param );
    };

  } // end namespace dgf

} // end namespace Dune

#endif
