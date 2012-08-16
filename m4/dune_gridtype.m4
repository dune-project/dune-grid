AC_DEFUN([DUNE_DEFINE_GRIDTYPE_INCLUDE],[dnl
m4_if($#,0,[],[dnl
  #include <$1>
m4_if($#,1,[],[DUNE_DEFINE_GRIDTYPE_INCLUDE(m4_shift($@))])dnl
])dnl
])


# DUNE_DEFINE_GRIDTYPE([GRIDTYPE],[ASSERTION],[DUNETYPE],[HEADER],...)
#
# Add a new GRIDTYPE target to DUNE's preprocessor magic.
# 
# Parameters: GRIDTYPE   name of the new target
#             ASSERTION  condition to be checked by the preprocessor
#             DUNETYPE   C++ type of the grid
#             HEADER     name of the header file which includes the grid
#
# Example: DUNE_DEFINE_GRIDTYPE([YASPGRID],[GRIDDIM == WORLDDIM],[Dune::YaspGrid< dimgrid >],[dune/grid/yaspgrid.hh],[dune/grid/io/file/dgfparser/dgfyasp.hh])
AC_DEFUN([DUNE_DEFINE_GRIDTYPE],[AH_BOTTOM(dnl
[/* add GRIDTYPE typedef for grid implementation $3:
    defining $1 during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if defined $1 && ! defined USED_$1_GRIDTYPE
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
DUNE_DEFINE_GRIDTYPE_INCLUDE(m4_shift(m4_shift(m4_shift($@))))dnl
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
