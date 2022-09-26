// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/yaspgrid.hh>

template<typename Grid>
bool test_gridinfo(const Grid& grid)
{
  Dune::gridinfo(grid);
  for (int level = 0; level <= grid.maxLevel(); ++level)
    Dune::gridlevellist(grid, level, "gridlevellist");
  Dune::gridleaflist(grid, "gridleaflist");

  return true;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  {
    Dune::YaspGrid<2> grid({1., 1.}, {4, 4});
    grid.globalRefine(2);
    test_gridinfo(grid);
  }

  {
    Dune::YaspGrid<3> grid({1., 1., 1.}, {4, 4, 4});
    grid.globalRefine(2);
    test_gridinfo(grid);
  }

  return 0;
}
