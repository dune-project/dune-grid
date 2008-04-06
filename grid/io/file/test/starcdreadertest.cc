// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <string>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/starcdreader.hh>
#include <dune/grid/test/gridcheck.cc>

template <class GridType>
void readGrid (const std::string& baseName)
{
  // read the grid
  std::auto_ptr<GridType> grid(Dune::StarCDReader<GridType>::read(baseName));

  std::cout << "Starting grid tests ." << std::flush;

  // check macro grid
  gridcheck(*grid);

  std::cout << " passed." << std::endl;
}


int main (int argc , char **argv) try
{
  std::string gridDirectory = "../../../../doc/grids/starcd/";

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
catch (...)
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
