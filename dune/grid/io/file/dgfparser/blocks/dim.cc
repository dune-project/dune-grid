// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/dim.hh>

namespace Dune
{

  namespace dgf
  {
    // DimBlock
    // --------

    DimBlock :: DimBlock ( std :: istream &in )
      : BasicBlock ( in, "Dimensions" )
    {
      if (isempty()) {
        DUNE_THROW(DGFException,
                   "no dimension of world specified!");
      } else {
        getnextline();
        line >> _dim;
        if (_dim<1) {
          DUNE_THROW(DGFException,
                     "negative dimension of world specified!");
        }
        else {
          if (noflines()==1)
            _dimworld=_dim;
          else {
            getnextline();
            line >> _dimworld;
            if (_dimworld < _dim) {
              DUNE_THROW(DGFException,
                         "negative dimension of world smaller than dim!");
            }
          }
        }
      }
    }

  } // end namespace dgf

} // end namespace Dune
