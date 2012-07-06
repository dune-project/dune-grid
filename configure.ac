#! /bin/bash
# $Id$
# Process this file with autoconf to produce a configure script.
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/grid/common/grid.hh])
AM_CONFIG_HEADER([config.h])

# --enable-extra-utils
DUNE_EXTRA_UTILS

# check all dune-module stuff
DUNE_CHECK_ALL

# set up flags for the automated test system
DUNE_AUTOBUILD_FLAGS

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS)')
AC_SUBST([LDADD], '$(top_builddir)/dune/grid/libgrid.la $(DUNE_LIBS)')

# Fix Albertas idea of where the dunegrid and libalbertagrid_?d libraries are,
# since it expects them in $(DUNE_GRID_LIBDIR) by default
AS_IF([test x"$HAVE_ALBERTA" = x"1"],
[
  ALBERTA_LIBS='$(top_builddir)/lib/libdunealbertagrid_$(ALBERTA_DIM)d.la $(top_builddir)/lib/libdunegrid.la $(ALBERTA_LIBPATHFLAGS) -lalberta_$(ALBERTA_DIM)d $(ALBERTA_BASE_LIBS)'
  ALBERTA1D_LIBS='$(top_builddir)/lib/libdunealbertagrid_1d.la $(top_builddir)/lib/libdunegrid.la $(ALBERTA_LIBPATHFLAGS) -lalberta_1d $(ALBERTA_BASE_LIBS)'
  ALBERTA2D_LIBS='$(top_builddir)/lib/libdunealbertagrid_2d.la $(top_builddir)/lib/libdunegrid.la $(ALBERTA_LIBPATHFLAGS) -lalberta_2d $(ALBERTA_BASE_LIBS)'
  ALBERTA3D_LIBS='$(top_builddir)/lib/libdunealbertagrid_3d.la $(top_builddir)/lib/libdunegrid.la $(ALBERTA_LIBPATHFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)'
])

DUNE_SUMMARY_ALL

# write output
AC_CONFIG_FILES([Makefile
    lib/Makefile
    doc/refinement/Makefile
    doc/doxygen/Makefile
    doc/doxygen/Doxyfile
    doc/grids/Makefile
    doc/grids/amc/Makefile
    doc/grids/amiramesh/Makefile
    doc/grids/dgf/Makefile
    doc/grids/gridfactory/Makefile
    doc/grids/gmsh/Makefile
    doc/grids/starcd/Makefile
    doc/Makefile
    dune/Makefile
    dune/grid/Makefile
    dune/grid/test/Makefile
    dune/grid/common/test/Makefile
    dune/grid/common/refinement/Makefile
    dune/grid/common/Makefile
    dune/grid/uggrid/Makefile
    dune/grid/sgrid/Makefile
    dune/grid/albertagrid/Makefile
    dune/grid/albertagrid/test/Makefile
    dune/grid/yaspgrid/Makefile
    dune/grid/onedgrid/Makefile
    dune/grid/utility/Makefile
    dune/grid/utility/test/Makefile
    dune/grid/alugrid/common/Makefile
    dune/grid/alugrid/2d/Makefile
    dune/grid/alugrid/3d/Makefile
    dune/grid/alugrid/Makefile
    dune/grid/io/file/test/Makefile
    dune/grid/io/file/amiramesh/Makefile
    dune/grid/io/file/dgfparser/Makefile
    dune/grid/io/file/dgfparser/blocks/Makefile
    dune/grid/io/file/dgfparser/test/Makefile
    dune/grid/io/file/dgfparser/utils/Makefile
    dune/grid/io/file/gnuplot/Makefile
    dune/grid/io/file/vtk/Makefile
    dune/grid/io/file/Makefile
    dune/grid/io/visual/grape/Makefile
    dune/grid/io/visual/Makefile
    dune/grid/io/Makefile
    dune/grid/geometrygrid/Makefile
    m4/Makefile
    dune-grid.pc
    src/Makefile
    src/alu-to-vtk/Makefile
    src/gmsh-to-alu/Makefile
    src/gridinfo-gmsh/Makefile
    src/structured-to-alu/Makefile
    ])
AC_OUTPUT
