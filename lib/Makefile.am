# $Id: $

if ALBERTA
ALBERTAGRID_LIBS = libdunealbertagrid_1d.la libdunealbertagrid_2d.la libdunealbertagrid_3d.la
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
libdunealbertagrid_1d_la_LIBADD =			\
	../dune/grid/albertagrid/libalbertagrid_1d.la	\
	libdunegrid.la

libdunealbertagrid_2d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_2d_la_SOURCES = dummy.cc
libdunealbertagrid_2d_la_LIBADD =			\
	../dune/grid/albertagrid/libalbertagrid_2d.la	\
	libdunegrid.la

libdunealbertagrid_3d_la_SOURCES =
nodist_EXTRA_libdunealbertagrid_3d_la_SOURCES = dummy.cc
libdunealbertagrid_3d_la_LIBADD =			\
	../dune/grid/albertagrid/libalbertagrid_3d.la	\
	libdunegrid.la

include $(top_srcdir)/am/global-rules
