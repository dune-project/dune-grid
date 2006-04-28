// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDTYPE_HH
#define DUNE_GRIDTYPE_HH
/**
 * @file
 * @brief  A simple strategy for defining a grid type depending on
 * defines set during the make proecess:
 *
 * @code
 *#if defined ALBERTAGRID && HAVE_ALBERTA
   typedef Dune::AlbertaGrid<dimworld,dimworld> GridType;
 *#elif defined ALUGRID_CUBE && HAVE_ALUGRID
   typedef Dune::ALU3dGrid<dimworld,dimworld,Dune::hexa> GridType;
 *#elif defined ALUGRID_SIMPLEX && HAVE_ALUGRID
   typedef Dune::ALU3dGrid<dimworld,dimworld,Dune::tetra> GridType;
 *#elif defined ONEDGRID
   typedef Dune::OneDGrid<1,1> GridType;
 *#elif defined SGRID
   typedef Dune::SGrid<dimworld,dimworld> GridType;
 *#elif defined YASPGRID
   typedef Dune::YaspGrid<dimworld,dimworld> GridType;
 *#else
 *#define YASPGRID
   typedef Dune::YaspGrid<dimworld,dimworld> GridType;
 *#endif
 * @endcode
 * The variable dimworld is determined by
 * @code
 *#if !defined GRIDDIM
   const int dimworld = DUNE_PROBLEM_DIM;
 *#else
 *#if !GRIDDIM
     const int dimworld = DUNE_PROBLEM_DIM;
 *#else
      const int dimworld = GRIDDIM;
 *#endif
 *#endif
 * @endcode
 * Remark:
 * -# By defauly Dune::YaspGrid<2,2> is used.
 * -# For \c ALBERTAGRID or \c ALUGRID with \c dimworld=1 Dune::OneDGrid is used.

 * To reduce differences between seriell and parallel runs as much as possible
 * additional macros are defined in this file
 * @code
 *#if HAVE_MPI_CPP &&
     (defined YASPGRID || defined  ALUGRID_SIMPLEX || defined ALUGRID_CUBE)
 *#include <mpi.h>
 *#define MPISTART \
   int myrank=-1; \
   int mysize = 1; \
   MPI_Init(&argc, &argv); \
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank); \
   MPI_Comm_size(MPI_COMM_WORLD,&mysize);
 *#define MPIEND \
   MPI_Finalize();
 *#else
 *#define MPISTART \
   int myrank=-1; \
   int mysize = 1;
 *#define MPI_COMM_WORLD -1
 *#define MPIEND
 *#endif
 * @endcode
 * To use this feature, MPISTART should always be called at the beginning
 * of the function main and MPIEND at the end of the main function:
 * @code
   int main(int argc, char ** argv, char ** envp) {
   MPISTART
   ...
   // construct the grid
   GridType* grid = MacroGrid(filename,MPI_COMM_WORLD);
   ...
   MPIEND
   }
 * @endcode
 * @author Andreas Dedner
 */


#if !defined GRIDDIM
const int dimworld = DUNE_PROBLEM_DIM;
  #define GRIDDIM DUNE_PROBLEM_DIM
  #warning --- No GRIDDIM defined, defaulting to DUNE_PROBLEM_DIM
#else
  #if !GRIDDIM
const int dimworld = DUNE_PROBLEM_DIM;
    #define GRIDDIM DUNE_PROBLEM_DIM
    #warning --- No GRIDDIM defined, defaulting to DUNE_PROBLEM_DIM
  #else
const int dimworld = GRIDDIM;
  #endif
#endif
#if defined ALBERTAGRID && HAVE_ALBERTA
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<1,1> GridType;
  #else
    #include "dgfalberta.hh"
typedef Dune::AlbertaGrid<dimworld,dimworld> GridType;
  #endif
const int refStepsForHalf = dimworld;
#elif defined ALUGRID_CUBE && HAVE_ALUGRID
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #include "dgfalu.hh"
typedef Dune::ALUCubeGrid<dimworld,dimworld> GridType;
  #endif
const int refStepsForHalf = 1;
#elif defined ALUGRID_SIMPLEX && HAVE_ALUGRID
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #include "dgfalu.hh"
typedef Dune::ALUSimplexGrid<dimworld,dimworld> GridType;
  #endif
const int refStepsForHalf = 1;
#elif defined ONEDGRID
  #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
const int refStepsForHalf = 1;
#elif defined SGRID
  #include "dgfs.hh"
typedef Dune::SGrid<dimworld,dimworld> GridType;
const int refStepsForHalf = 1;
#elif defined UGGRID
  #include "dgfug.hh"
typedef Dune::UGGrid<dimworld,dimworld> GridType;
const int refStepsForHalf = 1;
#elif defined YASPGRID
  #include "dgfyasp.hh"
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
const int refStepsForHalf = 1;
#else
  #define YASPGRID
  #include "dgfyasp.hh"
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
const int refStepsForHalf = 1;
  #warning --- No GRIDTYPE defined, defaulting to YASPGRID
#endif
#undef GRIDDIM

#if HAVE_MPI_CPP && (defined YASPGRID || defined  ALUGRID_SIMPLEX || defined ALUGRID_CUBE)
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

#endif
