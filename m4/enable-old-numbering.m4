AC_DEFUN([DUNE_OLD_NUMBERING],[
  # enable the old numbering methods in the grid interface
  AC_ARG_ENABLE(old-numbering,
   AC_HELP_STRING([--enable-old-numbering],[enable old numbering methods]),
   [enable_old_numbering=$enableval],
   [enable_old_numbering=yes])

  if test x$enable_old_numbering = xyes; then
    AC_DEFINE([DUNE_ENABLE_OLD_NUMBERING], [1], 
      [Define to 1 if the old numbering methods should be present.])
  fi
  AM_CONDITIONAL([DUNE_ENABLE_OLD_NUMBERING],
    [test x$enable_old_numbering = xyes])
])