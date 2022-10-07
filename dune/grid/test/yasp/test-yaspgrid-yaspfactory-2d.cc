// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/tensorgridfactory.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    std::string testID =
      "yaspfactory-2d-np" + std::to_string(mpiHelper.size());

    check_yasp(testID + "equidistant",
               YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true, 0));
    check_yasp(testID + "equidistantoffset",
               YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true, 0));
    check_yasp(testID + "tensor",
               YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(true, 0));

    // In 2D, also test refinement
    for (int refineOpt = 0; refineOpt <= 1; ++refineOpt) {
      std::string refTestID = testID + "-ref" + std::to_string(refineOpt);

      check_yasp(refTestID + "equidistant",
                 YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
      check_yasp(refTestID + "equidistantoffset",
                 YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
      check_yasp(refTestID + "tensor",
                 YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
    }

    // Test the generic constructor
    check_yasp(testID + "equidistant-generic-constructor",
               YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true,0,false,true));
    check_yasp(testID + "equidistantoffset-generic-constructor",
               YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true,0,false,true));
    check_yasp(testID + "tensor-generic-constructor",
               YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(true,0,false,true));

    // And periodicity
//    check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true, 0, true));
//    check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true, 0, true));
//    check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(true, 0, true));


  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
