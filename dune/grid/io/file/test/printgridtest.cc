// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/grid/yaspgrid.hh>

#include <dune/grid/io/file/printgrid.hh>

int main(int argc, char **argv)
{
  try {
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    // make grid
    const int dim = 2;
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim> N(Dune::fill_array<int,dim>(4));
    std::bitset<dim> periodic (false);
    periodic[0] = true;
    int overlap = 1;
    Dune::YaspGrid<dim> grid(L,N,periodic, overlap);

    // write .plt files (one for png, one for svg) without executing gnuplot on them
    Dune::printgrid (grid, "printgridtest_yasp_svg", 4000, false, false);
    Dune::printgrid (grid, "printgridtest_yasp_png", 4000, false);
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
