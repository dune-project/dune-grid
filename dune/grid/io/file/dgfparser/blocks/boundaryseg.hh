// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_BOUNDARYSEGBLOCK_HH
#define DUNE_DGF_BOUNDARYSEGBLOCK_HH

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <dune/grid/io/file/dgfparser/parser.hh>
#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {
    class BoundarySegBlock
      : public BasicBlock
    {
      int dimworld;                    // the dimension of the vertices (is given  from user)
      bool goodline;                   // active line describes a vertex
      std :: vector< unsigned int > p; // active vertex
      int bndid;
      typedef DGFBoundaryParameter::type BoundaryParameter;
      BoundaryParameter parameter;
      bool simplexgrid;

    public:
      typedef DGFEntityKey< unsigned int> EntityKey;
      typedef std::pair < int, BoundaryParameter > BndParam;

      // initialize vertex block and get first vertex
      BoundarySegBlock ( std :: istream &in, int pnofvtx,
                         int pdimworld, bool psimplexgrid );

      // some information
      int get( std :: map< EntityKey, BndParam > & facemap,
               bool fixedsize,
               int vtxoffset
               );

      bool ok()
      {
        return goodline;
      }

      int nofbound()
      {
        return noflines();
      }

    private:
      bool next();

      // get coordinates of active vertex
      int operator[] (int i)
      {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i<dimworld+1);
        return p[i];
      }

      int size()
      {
        return p.size();
      }

    };

  } // end namespace dgf

} // end namespace Dune

#endif
