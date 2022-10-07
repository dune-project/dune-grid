// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief A unit test for the TensorGridFactory
 */

#include <config.h>

#include <iostream>
#include <cassert>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/onedgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/utility/tensorgridfactory.hh>
#include <dune/grid/test/gridcheck.hh>

using namespace Dune;

template<class Grid>
void fillFactory(TensorGridFactory<Grid>& f)
{
  for (int i = 0; i<Grid::dimension; ++i)
  {
    f.setStart(i,-100.);
    f.fillIntervals(i,10,20.);
    f.fillRange(i, 5, 130.);
    f.geometricFillIntervals(i, 5, 2.0);
  }
}

int main (int argc , char **argv)
try {
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  // Test OneDGrid
  TensorGridFactory<OneDGrid> fac1;
  fillFactory(fac1);
  {
    auto grid = fac1.createGrid();
    gridcheck(*grid);
    std::cout << "OneDgrid with " << grid->size(0) << " cells created." << std::endl;
  }

#if HAVE_DUNE_UGGRID
  // Test UGGrid
  TensorGridFactory<UGGrid<2> > fac3;
  fillFactory(fac3);
  {
    auto grid = fac3.createGrid();
    gridcheck(*grid);
    std::cout << "UGGrid<2> with " << grid->size(0) << " cells created." << std::endl;
  }

  TensorGridFactory<UGGrid<3> > fac4;
  fillFactory(fac4);
  {
    auto grid = fac4.createGrid();
    gridcheck(*grid);
    std::cout << "UGGrid<3> with " << grid->size(0) << " cells created." << std::endl;
  }
#endif

  // Test YaspGrid
  TensorGridFactory<YaspGrid<2, TensorProductCoordinates<double, 2> > > fac5;
  fillFactory(fac5);
  {
    auto grid = fac5.createGrid();
    gridcheck(*grid);
    std::cout << "YaspGrid with " << grid->size(0) << " cells created." << std::endl;
  }

  return 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
