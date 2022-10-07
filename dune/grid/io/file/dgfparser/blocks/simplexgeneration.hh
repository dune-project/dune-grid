// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_SIMPLEXGENERATIONBLOCK_HH
#define DUNE_DGF_SIMPLEXGENERATIONBLOCK_HH

#include <iostream>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {

    class SimplexGenerationBlock
      : public BasicBlock
    {
      double area_;
      double angle_;
      bool display_;
      std::string path_;
      bool haspath_;
      std::string filename_;
      std::string filetype_;
      std::string parameter_;
      std::string dumpfilename_;
      bool hasfile_;
      int dimension_;

    public:
      SimplexGenerationBlock ( std :: istream &in );

      double maxArea ()
      {
        return area_;
      }

      double minAngle ()
      {
        return angle_;
      }

      bool display ()
      {
        return display_;
      }

      bool haspath ()
      {
        return haspath_;
      }

      std :: string path ()
      {
        return path_;
      }

      bool hasfile ()
      {
        return hasfile_;
      }

      std :: string filename ()
      {
        return filename_;
      }

      std :: string filetype ()
      {
        return filetype_;
      }

      int dimension ()
      {
        return dimension_;
      }

      std :: string parameter ()
      {
        return parameter_;
      }

      const std::string dumpFileName ( ) const
      {
        return dumpfilename_;
      }
    };

  } // end namespace dgf

} // end namespace Dune

#endif
