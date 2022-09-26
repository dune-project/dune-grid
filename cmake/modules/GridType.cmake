# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# Implement the GRIDTYPE preprocessor magic
#
# .. cmake_function:: dune_define_gridtype
#
#    .. cmake_param:: output
#       :single:
#       :required:
#       :positional:
#
#       String variable grid definition is written to.
#
#    .. cmake_param:: GRIDTYPE
#       :single:
#       :required:
#
#       The name of the grid type to register, .e.g. YASPGRID.
#
#    .. cmake_param:: DUNETYPE
#       :single:
#       :required:
#
#       The C++ type of the grid to be used for the typedef, e.g. YaspGrid< dimgrid >.
#
#    .. cmake_param:: ASSERTION
#       :single:
#
#       Condition to be checked by the preprocessor, e.g. GRIDDIM == WORLDDIM for grids like YaspGrid.
#
#    .. cmake_param:: HEADERS
#       :multi:
#       :required:
#
#       The header files that need to be included when using this grid.
#
#    This function registers a new type for the GRIDTYPE magic.
#

# option to enable GridSelector
option(DUNE_GRID_GRIDTYPE_SELECTOR "Grid selector definition added to config.h" OFF)

macro(dune_define_gridtype output)

cmake_parse_arguments(GRIDTYPE "" "GRIDTYPE;DUNETYPE;ASSERTION" "HEADERS" ${ARGN})
if(DUNE_GRID_GRIDTYPE_SELECTOR)
  if(NOT(GRIDTYPE_GRIDTYPE AND GRIDTYPE_DUNETYPE))
    message(ERROR "Both GRIDTYPE and DUNETYPE have to be set")
  endif(NOT(GRIDTYPE_GRIDTYPE AND GRIDTYPE_DUNETYPE))
  set(${output}
    "${${output}}
/* add GRIDTYPE typedef for grid implementation ${GRIDTYPE_DUNETYPE}:
   defining ${GRIDTYPE_GRIDTYPE} during compilation typedefs this grid implementation as GridType
   in namespace Dune::GridSelector;
   also integer constants dimgrid and dimworld are set in this namespace.
   The required headers for this grid implementation are also included.
*/
#if HAVE_DUNE_GRID && defined ${GRIDTYPE_GRIDTYPE} && ! defined USED_${GRIDTYPE_GRIDTYPE}_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error \"Ambiguous definition of GRIDTYPE.\"
  #endif

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error \"WORLDDIM < GRIDDIM does not make sense.\"
  #endif")
  if(GRIDTYPE_ASSERTION)
    set(${output} "${${output}}
  #if ! (${GRIDTYPE_ASSERTION})
     #error \"Preprocessor assertion ${GRIDTYPE_ASSERTION} failed.\"
  #endif"
      )
  endif(GRIDTYPE_ASSERTION)

  foreach(include ${GRIDTYPE_HEADERS})
    set(${output} "${${output}}
  #include <${include}>")
  endforeach(include ${GRIDTYPE_HEADERS})

  set(${output} "${${output}}
  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef ${GRIDTYPE_DUNETYPE} GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_${GRIDTYPE_GRIDTYPE}_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ${GRIDTYPE_GRIDTYPE} && ..")
# if disabled print message how to enable
else(DUNE_GRID_GRIDTYPE_SELECTOR)
  set(${output}
    "${${output}}\n/* ${GRIDTYPE_GRIDTYPE} not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */")
endif(DUNE_GRID_GRIDTYPE_SELECTOR)
endmacro(dune_define_gridtype output)
