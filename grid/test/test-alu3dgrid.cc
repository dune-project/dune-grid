// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include "gridcheck.cc"

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"

using namespace Dune;

template <class GridType>
void makeNonConfGrid(GridType &grid,int level,int adapt) {
  int myrank = grid.comm().rank();
  grid.loadBalance();
  grid.globalRefine(level);
  grid.loadBalance();
  for (int i=0; i<adapt; i++)
  {
    if (myrank==0)
    {
      typedef typename GridType :: template Codim<0> ::
      template Partition<Interior_Partition> :: LeafIterator LeafIterator;

      LeafIterator endit = grid.template leafend<0,Interior_Partition>   ();
      int nr = 0;
      int size = grid.size(0);
      for(LeafIterator it    = grid.template leafbegin<0,Interior_Partition> ();
          it != endit ; ++it,nr++ )
      {
        grid.mark(1,it);
        if (nr>size*0.8) break;
      }
    }
    grid.adapt();
    grid.loadBalance();
  }
}


template <class GridType>
void checkALUSerial(GridType & grid, int mxl = 2)
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
#if HAVE_MPI
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
  makeNonConfGrid(grid,gref,mxl);
  //gridcheck(grid);

  // -1 stands for leaf check
  checkCommunication(grid,grid.leafIndexSet(), -1, Dune::dvverb);

  for(int l=0; l<= mxl; ++l)
    checkCommunication(grid,grid.levelIndexSet(l), l , Dune::dvverb);
}
#else
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
  // only refine once
  checkALUSerial(grid,1);
}
#endif

int main (int argc , char **argv) {

  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpihelper = MPIHelper::instance(argc,argv);
  int myrank = mpihelper.rank();
  int mysize = mpihelper.size();

  try {
    /* use grid-file appropriate for dimensions */

    // extra-environment to check destruction
    {
      factorEpsilon = 5.e+5;
      // check empty grid

      if (myrank == 0)
        std::cout << "Check empty grids" << std::endl;
      {
        ALUCubeGrid<3,3> grid;
        checkALUSerial(grid);
      }
      {
        ALUSimplexGrid<3,3> grid;
        checkALUSerial(grid);
      }

      /*
         {
         std::string filename("alu-testgrid.triang");
         ALUSimplexGrid<2,2> grid(filename);
         checkALUSerial(grid,2);
         }
       */

      {
        std::string filename;
        if (mysize<=2)
          filename += "alu-testgrid.hexa";
        else
          filename += "largegrid_alu.hexa";

        ALUCubeGrid<3,3> grid (filename);

        if (myrank == 0) std::cout << "Check conform grid" << std::endl;
        checkALUParallel(grid,1,0);
        if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
        checkALUParallel(grid,0,2);
      }
      {
        std::string filename;
        if (mysize<=2)
          filename += "alu-testgrid.tetra";
        else
          filename += "examplegrid9.dgf.ALUgrid";

        ALUSimplexGrid<3,3> grid(filename);

        if (myrank == 0) std::cout << "Check conform grid" << std::endl;
        checkALUParallel(grid,0,0);  //1,3
        if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
        checkALUParallel(grid,0,2);  //1,3
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
