// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFGRIDTYPE_HH
#define DUNE_DGFGRIDTYPE_HH

/**
 * @file
 * @author Robert Kloefkorn
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
 *#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

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
 * programs \c CPPFLAGS one can choose the grids dimension and type by
 * invoking the following make command
 * @code
 *  make GRIDDIM=3 GRIDTYPE=ALBERTAGRID myprogram
 * @endcode
 * Here the value \c 3 for grid dimension and \c ALBERTAGRID for grid type are passed as
 * pre-processoer variables to gridtype.hh and are then used to define the \c typedef \c
 * GridType.
 *
 */

#include <dune/grid/utility/gridtype.hh>

// Note that we do not provide a default here, since gridtype.hh already already
// does this and sets the corresponding preprocessor flag. Since we only add
// some include files, we also need not check whether multiple grids are defined.

#if defined ALBERTAGRID
  #if not HAVE_ALBERTA
    #error "ALBERTAGRID defined but no ALBERTA version found!"
  #endif
//#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
  #include <dune/grid/albertagrid/dgfparser.hh>
#endif

#if defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALUGRID_CONFORM
  #if not HAVE_ALUGRID
    #error "ALUGRID_{CUBE,SIMPLEX,CONFORM} defined but no ALUGRID version found!"
  #endif
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif

#if defined ONEDGRID
  #include <dune/grid/io/file/dgfparser/dgfoned.hh>
#endif

#if defined SGRID
  #include <dune/grid/io/file/dgfparser/dgfs.hh>
#endif

#if defined UGGRID
  #if not HAVE_UG
    #error "UGGRID defined but no UG version found!"
  #endif
  #include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

#if defined YASPGRID
  #include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#endif

#endif // not defined DUNE_DGFGRIDTYPE_HH
