// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/simplexgeneration.hh>

namespace Dune
{

  namespace dgf
  {

    // SimplexGenerationBlock
    // ----------------------

    SimplexGenerationBlock :: SimplexGenerationBlock ( std :: istream &in )
      : BasicBlock(in,"Simplexgenerator"),
        area_(-1),
        angle_(-1),
        display_(false),
        haspath_(false),
        filetype_(),
        parameter_(),
        dumpfilename_(),
        hasfile_(false),
        dimension_(-1)
    {
      double x;
      bool b;
      int i;
      std::string p;
      if (findtoken("max-area"))
        if (getnextentry(x))
          area_=x;
      if (findtoken("min-angle"))
        if (getnextentry(x))
          angle_=x;
      if (findtoken("display"))
        if (getnextentry(b))
          display_=b;
      if (findtoken("path"))
        if (getnextentry(p)) {
          path_=p;
          haspath_=true;
        }
      if (findtoken("file")) {
        if (getnextentry(p)) {
          filename_=p;
          hasfile_=true;
        }
        if (getnextentry(p)) {
          filetype_=p;
        }
        if (findtoken("dimension"))
          if (getnextentry(i)) {
            dimension_=i;
          }
        gettokenparam("parameter",parameter_);
      }
      if (findtoken("dumpfilename"))
        if (getnextentry(p)) {
          dumpfilename_=p;
        }
    }

  } // end namespace dgf

} // end namespace Dune
