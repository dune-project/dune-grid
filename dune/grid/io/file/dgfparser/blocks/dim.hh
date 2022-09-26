// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_DIMBLOCK_HH
#define DUNE_DGF_DIMBLOCK_HH

#include <iostream>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {
    class DimBlock : public BasicBlock {
      int _dimworld;     // dimension of world
      int _dim;          // dimension of grid
    public:
      const static char* ID;
      // initialize block and get dimension of world
      DimBlock ( std :: istream &in );
      // get dimension of world found in block
      int dim() {
        return _dim;
      }
      int dimworld() {
        return _dimworld;
      }
      // some information
      bool ok() {
        return true;
      }
    };

  } // end namespace dgf

} // end namespace Dune

#endif
