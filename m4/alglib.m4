# AlgLib provides arbitrary precision linear algebra and quadratures (and more).
# see http://www.alglib.net/
# Unfortunatly the available downloads are rather buggy, therefore we have
# to provide a fixed version of the files used here - obtainable 
# by contacting dedner|nolte@mathematik.uni-freiburg.de

AC_DEFUN([DUNE_PATH_ALGLIB],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([DUNE_PATH_GMP])

  HAVE_ALGLIB=no

  AC_ARG_WITH(alglib,
    AS_HELP_STRING([--with-alglib=PATH],[directory to AlgLib for DUNE]))

  ac_save_PKG_CONFIG_PATH="$PKG_CONFIG_PATH"
  AS_IF([test x$with_alglib != x],[PKG_CONFIG_PATH="$with_alglib:$PKG_CONFIG_PATH"])
  AS_IF([pkg-config --atleast-version=1.0 alglib4dune],[
    HAVE_ALGLIB="yes (Version `pkg-config --modversion alglib4dune`)"
    ALGLIB_CPPFLAGS="`pkg-config --cflags alglib4dune` -DENABLE_ALGLIB=1"
    ALGLIB_LIBS="`pkg-config --libs alglib4dune`"
  ])
  PKG_CONFIG_PATH="$ac_save_PKG_CONFIG_PATH"

  AS_IF([test "$HAVE_ALGLIB" != no],[
    AC_LANG_PUSH([C++])
    ac_save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $ALGLIB_CPPFLAGS"
    AC_CHECK_HEADER([alglib/amp.h],[],[HAVE_ALGLIB=no])
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_LANG_POP
  ])

  AS_IF([test "$HAVE_ALGLIB" != no],[
    AC_DEFINE([HAVE_ALGLIB],[ENABLE_ALGLIB],[Was AlgLib for DUNE found and ALGLIB_CPPFLAGS used?])
    AC_SUBST([ALGLIB_CPPFLAGS],[$ALGLIB_CPPFLAGS])
    AC_SUBST([ALGLIB_LIBS],[$ALGLIB_LIBS])
    DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $ALGLIB_CPPFLAGS"
    DUNE_PKG_LIBS="$DUNE_PKG_LIBS $ALGLIB_LIBS"
  ])

  AM_CONDITIONAL(ALGLIB,[test "$HAVE_ALGLIB" != no])
  DUNE_ADD_SUMMARY_ENTRY([AlgLib for DUNE],[$HAVE_ALGLIB])
])
