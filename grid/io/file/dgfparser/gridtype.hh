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
 * additional the MPIHelper class is used to toggle between serial and parallel
 * runs via specialization of this class for runs using MPI und such that
 * doesn't.
 * To use this feature, the following code should always be called at the beginning
 * of the function main:
 * @code
   int main(int argc, char ** argv, char ** envp) {

   // get reference to the singelton MPIHelper
   MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);

   // optional one can get rank and size from this little helper class
   int myrank = mpiHelper.rank();
   int mysize = mpiHelper.size();

   ...
   // construct the grid
   GridPtr<GridType> grid(argv[1], mpiHelper.getCommunicator() );

   GridType & grid = *gridptr;
   ...

   // as the MPIHelper is a singleton, on it's destruction the
   // MPI_Finalize() command is called.
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

#endif
