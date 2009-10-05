#! /bin/bash
# $Id$
# Process this file with autoconf to produce a configure script.
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([grid/common/grid.hh])
AM_CONFIG_HEADER([config.h])

# make dune -> . symlink
DUNE_SYMLINK

# check all dune-module stuff
DUNE_CHECK_ALL

# set up flags for the automated test system
DUNE_AUTOBUILD_FLAGS

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS) $(DUNE_LIBS)')
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/grid/libgrid.la')

DUNE_SUMMARY_ALL

echo
echo Note: Most of the libraries checked for above are only used for the self-test
echo of Dune.  The library itself will build and the headers will work even if
echo Albert, MPI, etc. cannot be found.  An exception to this are UG and AmiraMesh
echo which need to be found right now if you want to use them later.
echo

# write output
AC_CONFIG_FILES([Makefile
    lib/Makefile
    doc/appl/refelements/Makefile
    doc/appl/Makefile
    doc/refinement/Makefile
    doc/doxygen/Makefile
    doc/doxygen/Doxyfile
    doc/Makefile
    grid/Makefile
    grid/test/Makefile
    grid/common/test/Makefile
    grid/common/quadraturerules/Makefile
    grid/common/refinement/Makefile
    grid/common/Makefile
    grid/uggrid/Makefile
    grid/sgrid/Makefile
    grid/albertagrid/Makefile
    grid/albertagrid/test/Makefile
    grid/yaspgrid/Makefile
    grid/onedgrid/Makefile
    grid/utility/Makefile
    grid/alugrid/2d/Makefile
    grid/alugrid/3d/Makefile
    grid/alugrid/Makefile
    grid/io/file/test/Makefile
    grid/io/file/amiramesh/Makefile
    grid/io/file/dgfparser/Makefile
    grid/io/file/dgfparser/test/Makefile
    grid/io/file/dgfparser/utils/Makefile
    grid/io/file/gnuplot/Makefile
    grid/io/file/vtk/Makefile
    grid/io/file/Makefile
    grid/io/visual/grape/Makefile
    grid/io/visual/Makefile
    grid/io/Makefile
    grid/genericgeometry/Makefile
    grid/genericgeometry/test/Makefile
    grid/geometrygrid/Makefile
    m4/Makefile
    dune-grid.pc])
AC_OUTPUT
