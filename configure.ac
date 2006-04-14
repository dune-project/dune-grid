#! /bin/bash
# $Id$
# Process this file with autoconf to produce a configure script.
AC_INIT(dune, 0.0, dune@hal.iwr.uni-heidelberg.de)
AM_INIT_AUTOMAKE(dune, 0.0, dune@hal.iwr.uni-heidelberg.de)
AC_CONFIG_SRCDIR([grid/common/grid.hh])
AM_CONFIG_HEADER([config.h])

# check all dune-module stuff
DUNE_CHECK_ALL_M

if test x$HAVE_DUNECOMMON != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-common module. Maybe you have to supply your dune-common directory as --with-dunecommon=dir])
fi

# preset variable to path such that #include <dune/...> works
AC_SUBST([AM_CPPFLAGS], '$(DUNECOMMON_CPPFLAGS) $(DUNEGRID_CPPFLAGS)')

AC_SUBST([LIBDUNE], '$(DUNECOMMONROOT)/common/libcommon.la $(DUNEGRIDROOT)/grid/libgrid.la')

DUNE_SUMMARY_ALL

echo
echo Note: Most of the libraries checked for above are only used for the self-test
echo of Dune.  The library itself will build and the headers will work even if
echo Albert, MPI, etc. cannot be found.  An exception to this are UG and AmiraMesh
echo which need to be found right now if you want to use them later.
echo

# write output
AC_CONFIG_FILES([Makefile
    doc/appl/refelements/Makefile
    doc/appl/Makefile
    doc/refinement/Makefile
    doc/doxygen/Makefile
    doc/layout/Makefile
    doc/Makefile
    grid/test/Makefile
    grid/common/test/Makefile
    grid/common/refinement/Makefile
    grid/common/Makefile
    grid/uggrid/Makefile
    grid/sgrid/Makefile
    grid/albertagrid/Makefile
    grid/yaspgrid/Makefile
    grid/onedgrid/Makefile
    grid/utility/Makefile
    grid/alu3dgrid/Makefile
    grid/Makefile
    io/file/test/Makefile
    io/file/amiramesh/Makefile
    io/file/alberta/Makefile
    io/file/Makefile
    io/visual/grape/Makefile
    io/visual/Makefile
    io/Makefile
    quadrature/test/Makefile
    quadrature/Makefile
    dune.pc])
AC_OUTPUT
