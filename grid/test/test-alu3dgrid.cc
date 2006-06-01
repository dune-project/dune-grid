// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

// #define ALUGRID_SIMPLEX
// #define GRIDDIM 3
// #include <dune/grid/io/file/dgfparser/gridtype.hh>
// #undef GRIDDIM
// #undef ALUGRID_SIMPLEX

#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#if HAVE_MPI_CPP
#include <mpi.h>
#define MPISTART \
  int myrank=-1; \
  int mysize = 1; \
  MPI_Init(&argc, &argv); \
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank); \
  MPI_Comm_size(MPI_COMM_WORLD,&mysize);
#define MPIEND \
  MPI_Finalize();
#else
#define MPISTART \
  int myrank=-1; \
  int mysize = 1;
#define MPI_COMM_WORLD -1
#define MPIEND
#endif

//#include <dune/grid/alugrid.hh>

#define ALUGRID_TESTING
#include "gridcheck.cc"
#undef ALUGRID_TESTING

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"


using namespace Dune;

template <class GridType>
void makeNonConfGrid(GridType &grid,int level,int adapt) {
  int myrank = grid.comm().rank();
  int mysize = grid.comm().size();
  grid.loadBalance();
  grid.globalRefine(level);
  grid.loadBalance();
  for (int i=0; i<adapt; i++) {
    if (myrank==0) {
      typedef typename GridType :: template Codim<0> ::
      template Partition<Interior_Partition> :: LeafIterator LeafIterator;
      LeafIterator endit = grid.template leafend<0,Interior_Partition>   ();
      int nr = 0;
      int size = grid.size(0);
      for(LeafIterator it    = grid.template leafbegin<0,Interior_Partition> ();
          it != endit ; ++it,nr++ ) {
        grid.mark(1,it);
        if (nr>size*0.8)
          break;
      }
    }
    grid.adapt();
    grid.loadBalance();
  }
}


template <class GridType>
void checkALUSerial(GridType & grid, int mxl = 2)
{
  return;
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
#if HAVE_MPI_CPP
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
  makeNonConfGrid(grid,gref,mxl);
  int myrank = grid.comm().rank();
  int mysize = grid.comm().size();
  // gridcheck(grid);
  checkCommunication(grid,gref,mxl,Dune::dvverb);
}
#else
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{}
#endif

int main (int argc , char **argv) {
  MPISTART
  try {
    /* use grid-file appropriate for dimensions */

    // extra-environment to check destruction
    {
      factorEpsilon = 5.e+5;
      // check empty grid
      if (myrank == 0)
        std::cout << "Check empty grids" << std::endl;
      {
        std::string filename("");
        ALUCubeGrid<3,3> grid(filename,MPI_COMM_WORLD);
        checkALUSerial(grid);
      }
      {
        std::string filename("");
        ALUSimplexGrid<3,3>
        grid(filename,MPI_COMM_WORLD);
        checkALUSerial(grid);
      }
      /*
         {
         std::string filename("alu-testgrid.triang");
         ALUSimplexGrid<2,2> grid(filename);
         checkALU(grid,0);
         }
       */
      {
        //std::string filename("alu-testgrid.hexa");
        std::string filename("largegrid_alu.hexa");
        ALUCubeGrid<3,3>* grid=new ALUCubeGrid<3,3>(filename,MPI_COMM_WORLD);
        //std::string filename("alu-testgrid.dgf");
        //GridPtr<GridType> grid(filename.c_str(),MPI_COMM_WORLD);
        if (myrank == 0)
          std::cout << "Check conform grid" << std::endl;
        checkALUParallel(*grid,1,0);
        if (myrank == 0)
          std::cout << "Check non-conform grid" << std::endl;
        checkALUParallel(*grid,0,2);
      }
      {
        // std::string filename("alu-testgrid.tetra");
        std::string filename("examplegrid9.dgf.ALUgrid");
        ALUSimplexGrid<3,3>* grid=new ALUSimplexGrid<3,3>(filename,MPI_COMM_WORLD);
        //std::string filename("examplegrid9.dgf");
        //GridPtr<GridType> grid(filename.c_str(),MPI_COMM_WORLD);
        if (myrank == 0)
          std::cout << "Check conform grid" << std::endl;
        checkALUParallel(*grid,0,0);  //1,3
        if (myrank == 0)
          std::cout << "Check non-conform grid" << std::endl;
        checkALUParallel(*grid,0,2);  //1,3
      }
    };

  } catch (Dune::Exception &e) {
    #if HAVE_MPI_CPP
    MPI_Finalize();
    #endif
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    #if HAVE_MPI_CPP
    MPI_Finalize();
    #endif
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  #if HAVE_MPI_CPP
  MPI_Finalize();
  #endif

  return 0;
}
