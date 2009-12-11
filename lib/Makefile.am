# $Id: $

#the dune-grid library
lib_LTLIBRARIES = libdunegrid.la

# this is just a renaming of libgrid to prevent name clashes
libdunegrid_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libdunegrid_la_SOURCES = dummy.cc
#libdunegrid_la_LIBADD = ../dune/grid/libgrid.la

include $(top_srcdir)/am/global-rules
