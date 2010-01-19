## -*- autoconf -*-
# $Id: alberta.m4 5156 2008-04-14 09:28:06Z christi $
# searches for alberta-headers and libs

# Substitutes the following make variables:
#   ALBERTA_DIM
#     value of --with-alberta-dim, or taken from --world-dim, or 2 per default
#
#   ALBERTAROOT = /usr/local/alberta
#     Root dir of the alberta installation.  Set from --with-alberta=...
#
#   ALBERTA_BASE_LIBS = $(ALBERTA_LIBPATHFLAGS) -lalberta_util $ALBERTA_EXTRA
#     LIBS that are always required independent of dimension
#
#   DUNEALBERTA_LIBPATHFLAGS = -L$(top_builddir)/lib
#     Library path required for libalbertagrid.  The above value is apropriate
#     when building dune-grid itself.  Other modules will call
#     DUNE_GRID_CHECK_MODULE to check for dune-grid, so DUNE_GRID_CHECK_MODULE
#     will reset DUNEALBERTA_LIBPATHFLAGS to something more apropriate.
#
#   ALBERTA_LIBPATHFLAGS = -L$(ALBERTAROOT)/lib
#     Library path required for alberta
#
#   ALBERTA_LIBS = $(DUNEALBERTA_LIBPATHFLAGS) -ldunealbertagrid_$(ALBERTA_DIM)d \
#              $(ALBERTA_LIBPATHFLAGS) -lalberta_$(ALBERTA_DIM)d \
#              $(ALBERTA_BASE_LIBS)
#     All LIBS required for the configured dimension
#
#   ALBERTA1D_LIBS, ALBERTA2D_LIBS, ALBERTA3D_LIBS
#     Likewise, but with the given dimension hardcoded.
#
#   ALBERTA_INCLUDE_CPPFLAGS = -I$(ALBERTAROOT)/include/alberta
#     Include path required for Alberta.
#
#   ALBERTA_DIM_CPPFLAGS = -DALBERTA_DIM=$(ALBERTA_DIM)
#     Pass the configured dimension as a CPP define
#
#   ALBERTA_CPPFLAGS = $(ALBERTA_INCLUDE_CPPFLAGS) $(ALBERTA_DIM_CPPFLAGS) -DENABLE_ALBERTA
#     All CPPFLAGS required for the configured dimension
#
#   ALBERTA1D_CPPFLAGS, ALBERTA2D_CPPFLAGS, ALBERTA3D_CPPFLAGS
#     Likewise, but with the given dimension hardcoded.
#
#   ALBERTA_LDFLAGS =
#     Just for consistency.
#
#   ALBERTA1D_LDFLAGS, ALBERTA2D_LDFLAGS, ALBERTA3D_LDFLAGS
#     Just for consistency.
#
#   If you want to use the the configured dimension, you have to use
#   $(ALBERTA_LIBS), $(ALBERTA_CPPFLAGS) and $(ALBERTA_LDFLAGS).
#
#   If want to use a specific dimension, say 2, you have to use
#   $(ALBERTA2D_LIBS), $(ALBERTA2D_CPPFLAGS) and $(ALBERTA2D_LDFLAGS).
#
# Defines the folling CPP macro
#   ALBERTA_DIM
#     The Alberta dimension this binary will be linked with.
#   DUNE_ALBERTA_VERSION
#     Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0
#   HAVE_ALBERTA
#     This is only true if alberta-library was found by configure 
#     _and_ if the application uses the ALBERTA_CPPFLAGS
#
# Defines the following automake conditional
#    ALBERTA
AC_DEFUN([DUNE_PATH_ALBERTA],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_F77])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([DUNE_PATH_OPENGL])

  ALBERTA_DIM='$(WORLDDIM)'

  AC_ARG_WITH(alberta,
    AC_HELP_STRING([--with-alberta=PATH],[directory where ALBERTA (ALBERTA
    version 2.0 or higher) is installed.  You can pass additional required
    libraries in the ALBERTA_EXTRA environment variable (in a form suitable
    for $LIBS)]))

  # do not use alberta debug lib 
  with_alberta_debug=0

  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  # LIBS=""

  ## do nothing if no --with-alberta was supplied
  if test x$with_alberta != x && test x$with_alberta != xno ; then

    # initialize to some default value
    ALBERTAROOT="/usr/local/alberta"
    # is --with-alberta=PATH used?
    if test -d $with_alberta ; then
      ALBERTAROOT=`cd $with_alberta && pwd`
    elif test "x$with_alberta" != "xyes" ; then
      AC_MSG_WARN([ALBERTA directory '$with_alberta' does not exist])
    fi

    ALBERTA_VERSION="2.0"

    # set variables so that tests can use them
    ALBERTA_INCLUDE_CPPFLAGS="-I$ALBERTAROOT/include/alberta"
    ALBERTA_DIM_CPPFLAGS='-DALBERTA_DIM=$(ALBERTA_DIM)'
    ALBERTA_CPPFLAGS='$(ALBERTA_INCLUDE_CPPFLAGS) $(ALBERTA_DIM_CPPFLAGS) -DENABLE_ALBERTA'
    ALBERTA1D_CPPFLAGS='$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=1 -DENABLE_ALBERTA'
    ALBERTA2D_CPPFLAGS='$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA'
    ALBERTA3D_CPPFLAGS='$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA'

    # check for header
    CPPFLAGS="$ac_save_CPPFLAGS $ALBERTA_INCLUDE_CPPFLAGS -DDIM_OF_WORLD=3 -DEL_INDEX=0"
    AC_CHECK_HEADER([alberta.h], [HAVE_ALBERTA="1"],
      AC_MSG_WARN([alberta.h not found in $ALBERTA_INCLUDE_CPPFLAGS]))

    if test "x$HAVE_ALBERTA" = "x1" ; then
      AC_CHECK_MEMBER([struct el_info.wall_bound],[ALBERTA_VERSION="3.0"],[AC_MSG_WARN([version 3 not found, now looking for version 2])],[#include <alberta.h>])
    fi

    CPPFLAGS="$ac_save_CPPFLAGS $ALBERTA_INCLUDE_CPPFLAGS"

    # TODO: check if static flag exists 
    # link_static_flag defines the flag for the linker to link only static
    # didnt work, with $link_static_flag, so quick hack here

    # check for libalberta_util...
    ALBERTA_LIBPATHFLAGS='-L$(ALBERTAROOT)/lib'
    DUNEALBERTA_LIBPATHFLAGS='-L$(top_builddir)/lib'
    LDFLAGS="$LDFLAGS -L$ALBERTAROOT/lib"
    if test "x$HAVE_ALBERTA" = "x1" ; then
      AC_CHECK_LIB(alberta_util,[alberta_calloc],
        [LIBS="-lalberta_util $LIBS"],
        [HAVE_ALBERTA="0"
         AC_MSG_WARN(-lalberta_util not found!)])
    fi

    # check for ALBERTA grid library...
    if test "x$HAVE_ALBERTA" = "x1" ; then

      # we do not check libraries for ALBERTA 3.0 (linking would require libtool)
      if test "$ALBERTA_VERSION" == "3.0" ; then
        AC_MSG_WARN([ALBERTA $ALBERTA_VERSION found -- Skipping check for ALBERTA grid libraries])
      else
        AC_CHECK_LIB(alberta_1d,[mesh_traverse], [],
          [HAVE_ALBERTA="0"
           AC_MSG_WARN(-lalberta_1d not found!)])
        AC_CHECK_LIB(alberta_2d,[mesh_traverse], [],
          [HAVE_ALBERTA="0"
           AC_MSG_WARN(-lalberta_2d not found!)])
        AC_CHECK_LIB(alberta_3d,[mesh_traverse], [],
          [HAVE_ALBERTA="0"
           AC_MSG_WARN(-lalberta_3d not found!)])
      fi
      ALBERTA_BASE_LIBS="\$(ALBERTA_LIBPATHFLAGS) -lalberta_util $ALBERTA_EXTRA"
      # define varaible lib name depending on problem and world dim, to change
      # afterwards easily 
      ALBERTA_LIBS='$(DUNEALBERTA_LIBPATHFLAGS) -ldunealbertagrid_$(ALBERTA_DIM)d $(ALBERTA_LIBPATHFLAGS) -lalberta_$(ALBERTA_DIM)d $(ALBERTA_BASE_LIBS)'
      ALBERTA1D_LIBS='$(DUNEALBERTA_LIBPATHFLAGS) -ldunealbertagrid_1d $(ALBERTA_LIBPATHFLAGS) -lalberta_1d $(ALBERTA_BASE_LIBS)'
      ALBERTA2D_LIBS='$(DUNEALBERTA_LIBPATHFLAGS) -ldunealbertagrid_2d $(ALBERTA_LIBPATHFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)'
      ALBERTA3D_LIBS='$(DUNEALBERTA_LIBPATHFLAGS) -ldunealbertagrid_3d $(ALBERTA_LIBPATHFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)'
    fi

  fi  # end of alberta check (--without wasn't set)

  # clear all the LDFLAGS variables explicitly
  ALBERTA_LDFLAGS=
  ALBERTA1D_LDFLAGS=
  ALBERTA2D_LDFLAGS=
  ALBERTA3D_LDFLAGS=

  # survived all tests?
  if test "x$HAVE_ALBERTA" = "x1" ; then
    AC_DEFINE(HAVE_ALBERTA, ENABLE_ALBERTA,
      [This is only true if alberta-library was found by configure 
       _and_ if the application uses the ALBERTA_CPPFLAGS])

    if test "$ALBERTA_VERSION" = "2.0" ; then
      AC_DEFINE([DUNE_ALBERTA_VERSION], [0x200], [Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0])
    elif test "$ALBERTA_VERSION" = "3.0" ; then
      AC_DEFINE([DUNE_ALBERTA_VERSION], [0x300], [Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0])
    else
      AC_MSG_ERROR([Internal Inconsistency: Invalid Alberta version reported: $ALBERTA_VERSION.])
    fi

    # add to global list
    DUNE_ADD_ALL_PKG([Alberta], [$ALBERTA_CPPFLAGS], [$ALBERTA_LDFLAGS], [$ALBERTA_LIBS])

    # set variable for summary
    with_alberta="yes (Version $ALBERTA_VERSION)"
  else
    # clear all variables
    ALBERTA_DIM= 
    ALBERTAROOT= 
    ALBERTA_BASE_LIBS=
    DUNEALBERTA_LIBPATHFLAGS=
    ALBERTA_LIBPATHFLAGS=
    ALBERTA_LIBS=
    ALBERTA1D_LIBS=
    ALBERTA2D_LIBS=
    ALBERTA3D_LIBS=
    ALBERTA_INCLUDE_CPPFLAGS=
    ALBERTA_DIM_CPPFLAGS=
    ALBERTA_CPPFLAGS=
    ALBERTA1D_CPPFLAGS=
    ALBERTA2D_CPPFLAGS=
    ALBERTA3D_CPPFLAGS=
    ALBERTA_LDFLAGS=
    ALBERTA1D_LDFLAGS=
    ALBERTA2D_LDFLAGS=
    ALBERTA3D_LDFLAGS=

    # set variable for summary
    with_alberta="no"
  fi
    
  AC_SUBST([ALBERTA_DIM]) 
  AC_SUBST([ALBERTAROOT]) 
  AC_SUBST([ALBERTA_BASE_LIBS])
  AC_SUBST([DUNEALBERTA_LIBPATHFLAGS])
  AC_SUBST([ALBERTA_LIBPATHFLAGS])
  AC_SUBST([ALBERTA_LIBS])
  AC_SUBST([ALBERTA1D_LIBS])
  AC_SUBST([ALBERTA2D_LIBS])
  AC_SUBST([ALBERTA3D_LIBS])
  AC_SUBST([ALBERTA_INCLUDE_CPPFLAGS])
  AC_SUBST([ALBERTA_DIM_CPPFLAGS])
  AC_SUBST([ALBERTA_CPPFLAGS])
  AC_SUBST([ALBERTA1D_CPPFLAGS])
  AC_SUBST([ALBERTA2D_CPPFLAGS])
  AC_SUBST([ALBERTA3D_CPPFLAGS])
  AC_SUBST([ALBERTA_LDFLAGS])
  AC_SUBST([ALBERTA1D_LDFLAGS])
  AC_SUBST([ALBERTA2D_LDFLAGS])
  AC_SUBST([ALBERTA3D_LDFLAGS])

  # also tell automake
  AM_CONDITIONAL(ALBERTA, test x$HAVE_ALBERTA = x1)

  # reset old values
  LIBS="$ac_save_LIBS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"

  DUNE_ADD_SUMMARY_ENTRY([ALBERTA],[$with_alberta])
])
