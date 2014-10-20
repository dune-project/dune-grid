#
# Defined macros:
#
# dune_define_gridtype(<output> GRIDTYPE <gridtype> DUNETYPE <dunetype>
#   [ASSERTION <assertion>] HEADERS header1 [ header2 ...])
#
# Addes a new GRIDTYPE target to DUNE's preprocessor magic.
# The generated magic will be stored in variable output
#
# Parameters: gridtype   name of the new target
#             assertion  condition to be checked by the preprocessor
#             dunetype   C++ type of the grid
#             header1     name of the header file which includes the grid
#
macro(dune_define_gridtype output)
  cmake_parse_arguments(GRIDTYPE "" "GRIDTYPE;DUNETYPE;ASSERTION" "HEADERS" ${ARGN})
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
  #ifndef HAVE_GRIDTYPE
    #define HAVE_GRIDTYPE 0
  #endif
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
  #undef HAVE_GRIDTYPE
  #define HAVE_GRIDTYPE 1
  #define USED_${GRIDTYPE_GRIDTYPE}_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ${GRIDTYPE_GRIDTYPE} && ..")
endmacro(dune_define_gridtype output)

