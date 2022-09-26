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
#include <dune/grid/yaspgrid/backuprestore.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    const auto & mpiHelper = Dune::MPIHelper::instance(argc, argv);

    std::string testID =
      "backuprestore-equidistantoffset-np" + std::to_string(mpiHelper.size());

    // check the backup restore facility
    check_backuprestore(testID,
                        YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid());

    // Test again with refinement
    check_backuprestore(testID + "-ref",
                        YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true, 1));

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
