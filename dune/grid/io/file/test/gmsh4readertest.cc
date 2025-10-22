// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmsh/gmsh4reader.hh>

#include <dune/grid/test/gridcheck.hh>

#include <dune/grid/uggrid.hh>

using namespace Dune;


int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  const std::string path = DUNE_GRID_EXAMPLE_GRIDS_PATH "gmsh/";

  { // Read a 2d ASCII test grid
    auto grid = Impl::Gmsh::Gmsh4Reader<UGGrid<2> >::createGridFromFile(path + "hybrid-testgrid-2d-v4-ascii.msh");

    gridcheck(*grid);
  }

  { // Read a 2d binary test grid
    auto grid = Impl::Gmsh::Gmsh4Reader<UGGrid<2> >::createGridFromFile(path + "hybrid-testgrid-2d-v4-binary.msh");

    gridcheck(*grid);
  }

  { // Read a 3d ASCII test grid
    auto grid = Impl::Gmsh::Gmsh4Reader<UGGrid<3> >::createGridFromFile(path + "hybrid-testgrid-3d-v4-ascii.msh");

    gridcheck(*grid);
  }

  { // Read a 3d binary test grid
    auto grid = Impl::Gmsh::Gmsh4Reader<UGGrid<2> >::createGridFromFile(path + "hybrid-testgrid-3d-v4-ascii.msh");

    gridcheck(*grid);
  }

  return 0;
}
