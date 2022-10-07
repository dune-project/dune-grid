// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil -*-

// Based on: Boilerplate tutorial poisson_uniform.cc

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <array>
#include <iostream>

#include <dune/grid/yaspgrid.hh>

#include <dune/grid/io/file/printgrid.hh>

int main(int argc, char **argv)
{
  if (std::system("gnuplot --version") != 0) {
    std::cerr << "GNUplot was not found." << std::endl;
    std::exit(77);
  }

  try {
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

    // make grid
    const int dim = 2;
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim> N;
    std::fill(N.begin(), N.end(), 4);
    std::bitset<dim> periodic (false);
    periodic[0] = true;
    int overlap = 1;
    Dune::YaspGrid<dim> grid(L,N,periodic, overlap);

    // write .plt files (one for png, one for svg) without executing gnuplot on them
    Dune::printGrid (grid, helper, "printgridtest_yasp_svg", 4000, true, false);
    Dune::printGrid (grid, helper, "printgridtest_yasp_png", 4000, true);
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  // done
  return 0;
}
