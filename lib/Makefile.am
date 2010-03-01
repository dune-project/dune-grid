# $Id: $

ALBERTAGRID_LIBS =

if ALBERTA_1D
ALBERTAGRID_LIBS += libdunealbertagrid_1d.la
endif
if ALBERTA_2D
ALBERTAGRID_LIBS += libdunealbertagrid_2d.la
endif
if ALBERTA_3D
ALBERTAGRID_LIBS += libdunealbertagrid_3d.la
endif
if ALBERTA_4D
ALBERTAGRID_LIBS += libdunealbertagrid_4d.la
endif
if ALBERTA_5D
ALBERTAGRID_LIBS += libdunealbertagrid_5d.la
endif
if ALBERTA_6D
ALBERTAGRID_LIBS += libdunealbertagrid_6d.la
endif
if ALBERTA_7D
ALBERTAGRID_LIBS += libdunealbertagrid_7d.la
endif
if ALBERTA_8D
ALBERTAGRID_LIBS += libdunealbertagrid_8d.la
endif
if ALBERTA_9D
ALBERTAGRID_LIBS += libdunealbertagrid_9d.la
endif


#the dune-grid library
lib_LTLIBRARIES = libdunegrid.la $(ALBERTAGRID_LIBS)

# this is just a renaming of libgrid to prevent name clashes
libdunegrid_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libdunegrid_la_SOURCES = dummy.cc
libdunegrid_la_LIBADD = ../dune/grid/libgrid.la

# ../dune/grid/albertagrid/libalbertagrid_?d.la is a convenience library, so
# its complete contents will be copied into libdunealbertagrid_?d.la.
# libdunegrid.la however is not a convenience library, so only the dependency
# will be recorded.

libdunealbertagrid_1d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_1d_la_SOURCES = dummy.cc
libdunealbertagrid_1d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_1d.la	libdunegrid.la

libdunealbertagrid_2d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_2d_la_SOURCES = dummy.cc
libdunealbertagrid_2d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_2d.la	libdunegrid.la

libdunealbertagrid_3d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_3d_la_SOURCES = dummy.cc
libdunealbertagrid_3d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_3d.la	libdunegrid.la

libdunealbertagrid_4d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_4d_la_SOURCES = dummy.cc
libdunealbertagrid_4d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_4d.la	libdunegrid.la

libdunealbertagrid_5d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_5d_la_SOURCES = dummy.cc
libdunealbertagrid_5d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_5d.la	libdunegrid.la

libdunealbertagrid_6d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_6d_la_SOURCES = dummy.cc
libdunealbertagrid_6d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_6d.la	libdunegrid.la

libdunealbertagrid_7d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_7d_la_SOURCES = dummy.cc
libdunealbertagrid_7d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_7d.la	libdunegrid.la

libdunealbertagrid_8d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_8d_la_SOURCES = dummy.cc
libdunealbertagrid_8d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_8d.la	libdunegrid.la

libdunealbertagrid_9d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_9d_la_SOURCES = dummy.cc
libdunealbertagrid_9d_la_LIBADD = ../dune/grid/albertagrid/libalbertagrid_9d.la	libdunegrid.la

include $(top_srcdir)/am/global-rules
