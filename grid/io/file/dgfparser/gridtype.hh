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
   typedef Dune::OneDGrid<1,1> GridType;
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
  #warning --- No GRIDDIM defined, defaulting to 2
const int dimworld = 2;
  #define GRIDDIM 2
#else
const int dimworld = GRIDDIM;
#endif

namespace Dune {
  /*! @brief Some simple static information for a given GridType
   */
  template <class GridType>
  struct DGFGridInfo {
    //! number of globalRefine steps needed to refuce h by 0.5
    static int refineStepsForHalf();
    //! relation between volume of children to volume of father.
    //! If this is not a constant the return value is -1
    static double refineWeight();
  };
}

#if defined ALBERTAGRID && HAVE_ALBERTA
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #include "dgfalberta.hh"
typedef Dune::AlbertaGrid<dimworld,dimworld> GridType;
  #endif
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return dimworld;}
    static double refineWeight() {return 0.5;}
  };
}
#elif defined ALUGRID_CUBE && HAVE_ALUGRID
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #if GRIDDIM == 2
      #error "ALUGRID_CUBE not implemented for dimension two!"
    #else
      #include "dgfalu.hh"
typedef Dune::ALUCubeGrid<dimworld,dimworld> GridType;
    #endif
  #endif
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dimworld);}
  };
}
#elif defined ALUGRID_SIMPLEX && HAVE_ALUGRID
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #include "dgfalu.hh"
typedef Dune::ALUSimplexGrid<dimworld,dimworld> GridType;
  #endif
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dimworld);}
  };
}
#elif defined ALUGRID_CONFORM && HAVE_ALUGRID
  #if GRIDDIM == 1
    #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
  #else
    #include "dgfalu.hh"
typedef Dune::ALUConformGrid<dimworld,dimworld> GridType;
  #endif
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return dimworld;}
    static double refineWeight() {return 0.5;}
  };
}
#elif defined ONEDGRID
  #include "dgfoned.hh"
typedef Dune::OneDGrid<dimworld,dimworld> GridType;
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return 0.5;}
  };
}
#elif defined SGRID
  #include "dgfs.hh"
typedef Dune::SGrid<dimworld,dimworld> GridType;
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dimworld);}
  };
}
#elif defined UGGRID && HAVE_UG
  #include "dgfug.hh"
typedef Dune::UGGrid<dimworld,dimworld> GridType;
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return -1.;}
  };
}
#elif defined YASPGRID
  #include "dgfyasp.hh"
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dimworld);}
  };
}
#else
// warning if grid type defined but package not found
  #if defined ALBERTAGRID
    #warning "ALBERTAGRID defined but no ALBERTA version found!"
  #endif
  #if defined ALUGRID_SIMPLEX || defined ALUGRID_CUBE
    #warning "ALUGRID_{SIMPLEX|CUBE} defined but no ALUGRID version found!"
  #endif
  #if defined UGGRID
    #warning "UGGRID defined but no UG version found!"
  #endif

//#define YASPGRID
  #include "dgfyasp.hh"
typedef Dune::YaspGrid<dimworld,dimworld> GridType;
  #warning --- No GRIDTYPE defined, defaulting to YASPGRID
namespace Dune {
  template <>
  struct DGFGridInfo<GridType> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dimworld);}
  };
}
#endif
#undef GRIDDIM
#endif
