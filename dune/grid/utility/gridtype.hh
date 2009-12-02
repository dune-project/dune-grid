// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDTYPE_HH
#define DUNE_GRIDTYPE_HH

/**
 * @file
 * @author Andreas Dedner
 *
 * @brief A simple strategy for defining a grid type depending on
 *        defines set during the make process.
 *
 * Check whether one of \c ALBERTAGRID, \c ALUGRID_CUBE, \c ALUGRID_SIMPLEX,
 * \c ALUGRID_CONFORM, \c ONEDGRID, \c SGRID, \c UGGRID or \c YASPGRID is
 * defined, and define the typedef \ref GridType accordingly.  Dimension of
 * grid and world are taken from dimgrid and dimworld.  If none of the above
 * are defined, default to \ref YaspGrid and generate a warning.  If more than
 * one of the above are defined, generate an error.  If the selected grid does
 * not support that combination of dimgrid and dimworld, generate an error.
 *
 * To reduce differences between serial and parallel runs as much as possible,
 * the Dune::MPIHelper class is used to toggle these runs.
 * To use this feature, the following code should always be called at the beginning
 * of the function main:
 * @code
 *#include <dune/grid/utility/gridtype.hh>

   ...

   int main(int argc, char ** argv, char ** envp) {

   // get reference to the singelton MPIHelper
   MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);

   // optional one can get rank and size from this little helper class
   int myrank = mpiHelper.rank();
   int mysize = mpiHelper.size();

   ...
   // construct the grid, see documentation for constructors
   GridType grid;
   ...

   // as the MPIHelper is a singleton, on it's destruction the
   // MPI_Finalize() command is called.
   }
 * @endcode
 *
 * A program that wants to use \ref gridtype.hh should be compiled with
 * either the <tt>ALL_PKG_CPPFLAGS</tt> or the <tt>GRIDDIM_CPPFLAGS</tt>
 * included in its <tt>CPPFLAGS</tt>.  It is then possible to specify the
 * desired grid at make time using something like
 * @code
 *  make GRIDDIM=3 GRIDTYPE=ALBERTAGRID myprogram
 * @endcode
 * However, for this to work, it is neccessary to give the
 * <tt>--with-grid-dim=...</tt> parameter at configure time.  In addition to
 * changing the default for \ref GRIDDIM (and thus for \ref dimgrid), this
 * parameter also enables the infrastructure in the makefiles by including the
 * following in <tt>ALL_PKG_CPPFLAGS</tt> or the <tt>GRIDDIM_CPPFLAGS</tt>:
 * @code
 *  -DGRIDDIM=$(GRIDDIM) -DWORLDDIM=$(WORLDDIM) -D$(GRIDTYPE)
 * @endcode
 */

#include <dune/grid/utility/griddim.hh>

// Check for AlbertaGrid
#if defined ALBERTAGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALBERTA
    #error "ALBERTAGRID defined but no ALBERTA version found!"
  #endif
  #if (GRIDDIM < 1) || (GRIDDIM > 3)
    #error "ALBERTAGRID is only available for GRIDDIM=1, GRIDDIM=2 and GRIDDIM=3."
  #endif

  #include <dune/grid/albertagrid.hh>
typedef Dune::AlbertaGrid< dimgrid > GridType;
  #define HAVE_GRIDTYPE 1
#endif


// Check ALUGrid
#if defined ALUGRID_CUBE
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_CUBE defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM != 3) || (WORLDDIM != GRIDDIM)
    #error "ALUGRID_CUBE is only available for GRIDDIM=3 and WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUCubeGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#if defined ALUGRID_SIMPLEX
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_SIMPLEX defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM < 2) || (GRIDDIM > 3)
    #error "ALUGRID_SIMPLEX is only available for GRIDDIM=2 and GRIDDIM=3."
  #endif
  #if (WORLDDIM != GRIDDIM)
    #error "ALUGRID_SIMPLEX is only available for WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUSimplexGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#if defined ALUGRID_CONFORM
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_CONFORM defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM != 2) || (WORLDDIM != GRIDDIM)
    #error "ALUGRID_CONFORM is only available for GRIDDIM=2 and WORLDDIM=GRIDDIM."
  #endif
  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUConformGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif


// Check OneDGrid
#if defined ONEDGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if (GRIDDIM != 1) || (WORLDDIM != GRIDDIM)
    #error "ONEDGRID is only available for GRIDDIM=1 and WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/onedgrid.hh>
typedef Dune :: OneDGrid GridType;
  #define HAVE_GRIDTYPE 1
#endif


// Check SGrid
#if defined SGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif

  #include <dune/grid/sgrid.hh>
typedef Dune :: SGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif


// Check UGGrid
#if defined UGGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_UG
    #error "UGGRID defined but no UG version found!"
  #endif
  #if (GRIDDIM < 2) || (GRIDDIM > 3)
    #error "UGGRID is only available for GRIDDIM=2 and GRIDDIM=3."
  #endif
  #if (GRIDDIM != WORLDDIM)
    #error "UGGRID only supports GRIDDIM=WORLDDIM."
  #endif

  #include <dune/grid/uggrid.hh>
typedef Dune :: UGGrid< dimgrid > GridType;
  #define HAVE_GRIDTYPE 1
#endif


// Check YASPGrid
#if defined YASPGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if (GRIDDIM != WORLDDIM)
    #error "YASPGRID only supports GRIDDIM=WORLDDIM."
  #endif

  #include <dune/grid/yaspgrid.hh>
typedef Dune :: YaspGrid< dimgrid > GridType;
  #define HAVE_GRIDTYPE 1
#endif


// default grid type
#ifndef HAVE_GRIDTYPE
  #if (GRIDDIM != WORLDDIM)
    #error "No default grid available for GRIDDIM<WORLDDIM."
    #define HAVE_GRIDTYPE 0
  #else
    #warning "No GRIDTYPE defined, defaulting to YASPGRID."

//! \cond
    #define YASPGRID
//! \endcond
    #include <dune/grid/yaspgrid.hh>
//! Default grid type
/**
 * Set with <tt>configure --with-grid-type</tt>.  If unset, defaults to
 * YaspGrid<dimgrid>.
 *
 * \note This is in general not the type of "the" grid, it only provides a
 *       default for certain cases.
 */
typedef Dune :: YaspGrid< dimgrid > GridType;
//! \cond
    #define HAVE_GRIDTYPE 1
//! \endcond
  #endif
#endif

#endif
