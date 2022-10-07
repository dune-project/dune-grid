// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/tensorgridfactory.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);

    // check the factory class for tensorproduct grids
    Dune::TensorGridFactory<Dune::YaspGrid<2, Dune::TensorProductCoordinates<double,2> > > factory;
    factory.setStart(0,-100.);
    factory.fillIntervals(0,10,20.);
    factory.fillRange(0, 5, 130.);
    factory.geometricFillIntervals(0, 5, 2.0);

    factory.geometricFillRange(1,10,100.,1.,false);
    factory.fillRange(1,10,200);
    factory.geometricFillRange(1,10,250.,1.,true);
    factory.fillUntil(1,50,1000.);

    factory.createGrid();

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
