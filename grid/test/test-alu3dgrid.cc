// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/grid/alugrid.hh>

#define ALUGRID_TESTING
#include "gridcheck.cc"
#undef ALUGRID_TESTING

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"

using namespace Dune;


template <class GridType>
void checkALU(GridType & grid, int mxl = 2)
{
  // be careful, each global refine create 8 x maxlevel elements
  gridcheck(grid);
  for(int i=0; i<mxl; i++) {
    grid.globalRefine(1);
    gridcheck(grid);
  }

  // check the method geometryInFather()
  checkGeometryInFather(grid);

  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
}

int main () {
  try {
    /* use grid-file appropriate for dimensions */

    // extra-environment to check destruction
    {
      factorEpsilon = 500.0;

      {
        std::string filename("alu-testgrid.hexa");
        ALUCubeGrid<3,3> grid(filename);
        checkALU(grid);
      }

      {
        std::string filename("alu-testgrid.tetra");
        ALUSimplexGrid<3,3>
        grid(filename);
        checkALU(grid);
      }

      {
        std::string filename("alu-testgrid.triang");
        ALUSimplexGrid<2,2> grid(filename);
        checkALU(grid,0);
      }
    };

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
