# define GRIDDIM_CPPFLAGS and add to DUNE_PKG_CPPFLAGS
# This defines GRIDDIM, WORLDDIM, and GRIDTYPE and assigns invalid values 
AC_DEFUN([DUNE_GRID_DIMENSION],[
  griddim_cppflags="-DGRIDDIM=$``(``GRIDDIM``)`` -DWORLDDIM=$``(``WORLDDIM``)`` -D$``(``GRIDTYPE``)``"
  AC_SUBST(GRIDDIM, 0)
  AC_SUBST(WORLDDIM, "$``(``GRIDDIM``)``")

  AC_SUBST(GRIDTYPE, [NOGRID])
  AC_SUBST(GRIDDIM_CPPFLAGS, $griddim_cppflags)
  DUNE_ADD_ALL_PKG([GRIDDIM], [$griddim_cppflags])
  # AC_MSG_RESULT([yes (GRIDDIM=$GRIDDIM, WORLDDIM=GRIDDIM and GRIDTYPE=$GRIDTYPE)])
])
