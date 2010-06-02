AC_DEFUN([DUNE_DEFINE_GRIDTYPE_INCLUDE],[dnl
m4_if(regexp([$1],'\w'),[-1],[dnl
  #include <$1>
],[dnl
  #include <m4_substr([$1],0,regexp([$1],'\w'))>
  DUNE_DEFINE_GRIDTYPE_INCLUDE([m4_substr([$1],[regexp([$1],'\w')])])dnl
])dnl
])

# DUNE_DEFINE_GRIDTYPE([GRIDTYPE],[ASSERTION],[DUNETYPE],[HEADER],[DGFHEADER])
#
# Add a new GRIDTYPE target DUNE's preprocessor magic.
# 
# Parameters: GRIDTYPE   name of the new target
#             ASSERTION  condition to be checked by the preprocessor
#             DUNETYPE   C++ type of the grid
#             HEADER     name of the header file which includes the grid
#             DGFHEADER  name of the header for the DGFGridFactory for the grid
#
# Example: DUNE_DEFINE_GRIDTYPE([YASPGRID],[GRIDDIM == WORLDDIM],[Dune::YaspGrid< dimgrid >],[dune/grid/yaspgrid.hh],[dune/grid/io/file/dgfparser/dgfyasp.hh])
AC_DEFUN([DUNE_DEFINE_GRIDTYPE],[AH_BOTTOM(dnl
[#if defined $1 && ! defined USED_$1_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif
]dnl
m4_if([$2],[],[],[
  #if ! ($2)
    #error "Preprocessor assertion $2 failed."
  #endif
])
DUNE_DEFINE_GRIDTYPE_INCLUDE([$4])dnl
DUNE_DEFINE_GRIDTYPE_INCLUDE([$5])dnl
[
  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef $3 GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_$1_GRIDTYPE 1
#endif // #if defined $1]dnl
)])
