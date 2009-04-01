# $Id$

if ALBERTA
noinst_LTLIBRARIES = libalbertagrid_1.la libalbertagrid_2.la libalbertagrid_3.la

libalbertagrid_1_la_SOURCES = elementinfo.cc
libalbertagrid_1_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DALBERTA_DIM=1
libalbertagrid_1_la_LIBADD = $(DUNE_COMMON_LIBS)

libalbertagrid_2_la_SOURCES = elementinfo.cc
libalbertagrid_2_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DALBERTA_DIM=2
libalbertagrid_2_la_LIBADD = $(DUNE_COMMON_LIBS)

libalbertagrid_3_la_SOURCES = elementinfo.cc
libalbertagrid_3_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DALBERTA_DIM=3
libalbertagrid_3_la_LIBADD = $(DUNE_COMMON_LIBS)
endif

albertadir = $(includedir)/dune/grid/albertagrid/
alberta_HEADERS = agrid.hh albertagrid.cc \
  agmemory.hh albertaextra.hh albertaheader.hh \
  indexsets.hh datahandle.hh griddim.hh \
  misc.hh  elementinfo.hh  meshpointer.hh macrodata.hh \
  dofadmin.hh  dofvector.hh  refinement.hh  coordcache.hh  level.hh \
  undefine-2.0.hh  undefine-2.1.hh \
  entity.hh  entity.cc  entitypointer.hh  hierarchiciterator.hh \
  algebra.hh  geometry.hh  geometry.cc  projection.hh  transformation.hh \
  leveliterator.hh  leafiterator.hh  treeiterator.hh  treeiterator.cc \
  intersection.hh  intersection.cc  intersectioniterator.hh \
  capabilities.hh  gridfactory.hh  dgfparser.hh  albertareader.hh

SUBDIRS = test

include $(top_srcdir)/am/global-rules
