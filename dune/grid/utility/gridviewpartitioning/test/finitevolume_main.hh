// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_MAIN_HH
#define DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_MAIN_HH

#include <algorithm>
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <sstream>
#include <vector>                 // STL vector class

#include <mpi.h>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/utility/structuredgridfactory.hh>

#include "vtkout.hh"
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

template<class G, class Evolve, class EvolveOnIntriorIntersection>
void timeloop (const std::string &desc, const std::string &prefix,
               const G& grid, Evolve evolve, double tend,
               const EvolveOnIntriorIntersection &evolveOnII)
{
  std::cout << "FV (" << desc << "):" << std::endl;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,Dune::MCMGElementLayout>
  mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);                           /*@\label{fvc:init}@*/
  vtkout(grid, c, prefix.c_str(), 0, 0.0);

  // now do the time steps
  double t=0,dt;
  int k=0;
  // disable saving inside the loop.
  const double saveInterval = tend;
  double saveStep = saveInterval;
  int counter = 1;

  double start = MPI_Wtime();
  while (t<tend)                                       /*@\label{fvc:loop0}@*/
  {
    // augment time step counter
    ++k;

    // apply finite volume scheme
    evolve(grid,mapper,c,t,dt, evolveOnII);

    // augment time
    t += dt;

    // check if data should be written
    if (t >= saveStep && t < tend)
    {
      // write data
      vtkout(grid, c, prefix.c_str(), counter, t);

      // increase counter and saveStep for next interval
      saveStep += saveInterval;
      ++counter;
    }

    // print info about time, timestep size and counter
    // std::cout << "s=" << grid.size(0)
    //           << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
  }                                                    /*@\label{fvc:loop1}@*/
  double end = MPI_Wtime();
  std::cout << "Elapsed: " << (end - start) << " s" << std::endl;

  // output results
  vtkout(grid, c, prefix.c_str(), counter, tend);     /*@\label{fvc:file}@*/
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    Dune::FieldVector<Grid::ctype, Grid::dimensionworld> lowerLeft(0);
    Dune::FieldVector<Grid::ctype, Grid::dimensionworld> upperRight(1);
    Dune::array<unsigned,  Grid::dimension> elements;
    std::fill(elements.begin(), elements.end(), 1);
    Dune::shared_ptr<Grid> gridp;
    if(Dune::Capabilities::hasSingleGeometryType<Grid>::v &&
       Dune::GeometryType
         (Dune::Capabilities::hasSingleGeometryType<Grid>::topologyId,
          Grid::dimension).isSimplex())
    {
      gridp = Dune::StructuredGridFactory<Grid>::
        createSimplexGrid(lowerLeft, upperRight, elements);
    }
    else
    {
      gridp = Dune::StructuredGridFactory<Grid>::
        createCubeGrid(lowerLeft, upperRight, elements);
    }

    // grid reference
    Grid& grid = *gridp;

    int level = 0;
    if(argc >= 2) {
      std::istringstream s(argv[1]);
      s >> level;
      if(s.fail()) {
        std::cerr << "Invalid argument: " << argv[1] << std::endl;
        return 1;
      }
    }

    // refine grid until upper limit of level
    grid.globalRefine(level);

    std::cout << "Timer Ganularity: " << MPI_Wtick() << std::endl;
    std::cout << "Capabilities::viewThreadSafe<" << gridDescription << ">::v=="
              << Dune::Capabilities::viewThreadSafe<Grid>::v << std::endl;

    // do time loop until end time 0.5
    timeloop(gridDescription+", Seqential/Plain",
             gridPrefix+"_seq_plain_concentration", grid, SeqEvolve(), 0.5,
             EvolveOnInteriorIntersectionPlain());
    // do time loop until end time 0.5
    timeloop(gridDescription+", Seqential/Opt",
             gridPrefix+"_seq_opt_concentration", grid, SeqEvolve(), 0.5,
             EvolveOnInteriorIntersectionOptimized());
    // do time loop until end time 0.5
    timeloop(gridDescription+", OSThreads(2)/Plain",
             gridPrefix+"_ost_plain_concentration", grid, OSThreadsEvolve(), 0.5,
             EvolveOnInteriorIntersectionPlain());
    // do time loop until end time 0.5
    timeloop(gridDescription+", OSThreads(2)/Opt",
             gridPrefix+"_ost_opt_concentration", grid, OSThreadsEvolve(), 0.5,
             EvolveOnInteriorIntersectionOptimized());
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}

#endif // DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_MAIN_HH
