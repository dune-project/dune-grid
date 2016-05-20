// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <config.h>

#include <string>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/starcdreader.hh>
#include <dune/grid/test/gridcheck.hh>

template <class GridType>
void readGrid (const std::string& baseName)
{
  // read the grid
  std::unique_ptr<GridType> grid(Dune::StarCDReader<GridType>::read(baseName));

  std::cout << "Starting grid tests ." << std::flush;

  // check macro grid
  gridcheck(*grid);

  std::cout << " passed." << std::endl;
}


int main (int argc , char **argv) try
{
  // initialize MPI if necessary
  Dune::MPIHelper::instance(argc, argv);

  std::string gridDirectory = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "starcd/";

  // Currently, Star-CD only works for UGGrid
  readGrid<Dune::UGGrid<3> >(gridDirectory + "star");
  readGrid<Dune::UGGrid<3> >(gridDirectory + "tets");
  readGrid<Dune::UGGrid<3> >(gridDirectory + "withprism");
  readGrid<Dune::UGGrid<3> >(gridDirectory + "withpyramid");

  return 0;
}
catch (Dune::Exception& e)
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e)
{
  std::cerr << e.what() << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
