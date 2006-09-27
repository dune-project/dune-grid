# $Id$

libalberta_la_CXXFLAGS = @AM_CPPFLAGS@ $(ALBERTA_CPPFLAGS)
libalberta_la_LIBADD =

if ALBERTA
  noinst_LTLIBRARIES = libalberta.la
  libalberta_la_SOURCES = albertaextra.cc
#else
#  noinst_LTLIBRARIES = libalu3d.la
#  libalu3d_la_SOURCES = topology.cc
endif

albertadir = $(includedir)/dune/grid/albertagrid/
alberta_HEADERS = agrid.hh agelementindex.cc albertagrid.cc \
  agcommunicator.hh agmemory.hh albertaextra.hh albertaheader.hh \
  alberta_undefs.hh indexsets.hh referencetopo.hh datahandle.hh

include $(top_srcdir)/am/global-rules
