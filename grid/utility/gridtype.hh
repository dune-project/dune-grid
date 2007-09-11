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
   typedef Dune::ALUCubeGrid<dimworld,dimworld> GridType;
 *#elif defined ALUGRID_SIMPLEX && HAVE_ALUGRID
   typedef Dune::ALUSimplexGrid<dimworld,dimworld> GridType;
 *#elif defined ONEDGRID
   typedef Dune::OneDGrid GridType;
 *#elif defined SGRID
   typedef Dune::SGrid<dimworld,dimworld> GridType;
 *#elif defined YASPGRID
   typedef Dune::YaspGrid<dimworld,dimworld> GridType;
 *#else
   // default GridType is YaspGrid
   typedef Dune::YaspGrid<dimworld,dimworld> GridType;
 *#endif
 * @endcode
 * The variable dimworld is determined by
 * @code
   // default value is 2
   const int dimworld = GRIDDIM;
 * @endcode
 * Remark:
 * -# By defauly Dune::YaspGrid<2,2> is used.
 * -# For \c ALBERTAGRID or \c ALUGRID with \c dimworld=1 Dune::OneDGrid is used.

 * To reduce differences between seriell and parallel runs as much as possible
 * additional the Dune::MPIHelper class is used to toggle between serial and parallel
 * runs via specialization of this class for runs using MPI und such that
 * doesn't.
 * To use this feature, the following code should always be called at the beginning
 * of the function main:
 * @code
 *#include <dune/grid/io/file/dgfparser/gridtype.hh>

   ...

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
 * To make use of this feature, one has to configure DUNE by using the
 * \c --with-grid-dim=1|2|3 and optional --with-grid-type=ALBERTAGRID |
 * ALUGRID_CUBE | ALUGRID_SIMPLEX | ONEDGRID | SGRID | UGGRID | YASPGRID.
 * This will add the following to the ALL_PKG_CPPFLAGS
 * @code
 * -DGRIDDIM=$(GRIDDIM) -DGRIDTYPE=$(GRIDTYPE)
 * @endcode
 * No by adding the \c ALL_PKG_CPPFLAGS or \c GRIDDIM_CPPFLAGS to the
 * programs \c CXXFLAGS one can choose the grids dimension and type by
 * invoking the following make command
 * @code
 *  make GRIDDIM=3 GRIDTYPE=ALBERTAGRID myprogram
 * @endcode
 * Here the value \c 3 for grid dimension and \c ALBERTAGRID for grid type are passed as
 * pre-processoer variables to gridtype.hh and are then used to define the \c typedef \c
 * GridType.
 *
 * @author Andreas Dedner
 */

#ifndef GRIDDIM
  #warning --- No GRIDDIM defined, defaulting to 3
const int dimworld = 3;
  #define GRIDDIM 3
#else
const int dimworld = GRIDDIM;
#endif

#if defined ALBERTAGRID
  #if not HAVE_ALBERTA
    #error "ALBERTAGRID defined but no ALBERTA version found!"
  #endif
  #if GRIDDIM < 2 || GRIDDIM > 3
    #error "ALBERTAGRID is only available for GRIDDIM=2 and GRIDDIM=3"
  #endif
  #include <dune/grid/albertagrid.hh>
typedef Dune::AlbertaGrid<dimworld,dimworld> GridType;
#elif defined ALUGRID_CUBE
  #if not HAVE_ALUGRID
    #error "ALUGRID_CUBE defined but no ALUGRID version found!"
  #endif
  #if GRIDDIM != 3
    #error ALUGRID_CUBE is only available for GRIDDIM=3
  #endif
  #include <dune/grid/alugrid.hh>
typedef Dune::ALUCubeGrid<dimworld,dimworld> GridType;
#elif defined ALUGRID_SIMPLEX
  #if not HAVE_ALUGRID
    #error "ALUGRID_SIMPLEX defined but no ALUGRID version found!"
  #endif
  #if GRIDDIM < 2 || GRIDDIM > 3
    #error ALUGRID_SIMPLEX is only available for GRIDDIM=2 and GRIDDIM=3
  #endif
  #include <dune/grid/alugrid.hh>
typedef Dune::ALUSimplexGrid<dimworld,dimworld> GridType;
#elif defined ALUGRID_CONFORM && HAVE_ALUGRID
  #if not HAVE_ALUGRID
    #error "ALUGRID_CONFORM defined but no ALUGRID version found!"
  #endif
  #if GRIDDIM != 3
    #error ALUGRID_CONFORM is only available for GRIDDIM=3
  #endif
  #include <dune/grid/alugrid.hh>
typedef Dune::ALUConformGrid<dimworld,dimworld> GridType;
#elif defined ONEDGRID
  #if GRIDDIM != 1
    #error ONEDGRID is only available for GRIDDIM=1
  #endif
  #include <dune/grid/onedgrid.hh>
typedef Dune::OneDGrid GridType;
#elif defined SGRID
  #include <dune/grid/sgrid.hh>
typedef Dune::SGrid<dimworld,dimworld> GridType;
#elif defined UGGRID
  #if not HAVE_UG
    #error "UGGRID defined but no UG version found!"
  #endif
  #if GRIDDIM < 2 || GRIDDIM > 3
    #error UGGRID is only available for GRIDDIM=2 and GRIDDIM=3
  #endif
  #include <dune/grid/uggrid.hh>
typedef Dune::UGGrid<dimworld> GridType;
#elif defined YASPGRID
  #include <dune/grid/yaspgrid.hh>
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
#else
// fallback
  #include <dune/grid/yaspgrid.hh>
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
  #warning --- No GRIDTYPE defined, defaulting to YASPGRID
#endif
#undef GRIDDIM
#endif
