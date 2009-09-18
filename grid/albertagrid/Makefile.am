# $Id$

if ALBERTA
noinst_LTLIBRARIES = libalbertagrid_1.la libalbertagrid_2.la libalbertagrid_3.la

libalbertagrid_1_la_SOURCES = elementinfo.cc  meshpointer.cc
libalbertagrid_1_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=1 -DENABLE_ALBERTA
libalbertagrid_1_la_LIBADD = $(DUNE_COMMON_LIBS)

libalbertagrid_2_la_SOURCES = elementinfo.cc  meshpointer.cc
libalbertagrid_2_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
libalbertagrid_2_la_LIBADD = $(DUNE_COMMON_LIBS)

libalbertagrid_3_la_SOURCES = elementinfo.cc  meshpointer.cc
libalbertagrid_3_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
libalbertagrid_3_la_LIBADD = $(DUNE_COMMON_LIBS)
endif

albertadir = $(includedir)/dune/grid/albertagrid/
alberta_HEADERS = agrid.hh albertagrid.cc \
  agmemory.hh  albertaheader.hh \
  indexsets.hh datahandle.hh \
  misc.hh  macroelement.hh  elementinfo.hh  meshpointer.hh macrodata.hh \
  dofadmin.hh  dofvector.hh  refinement.hh  coordcache.hh  level.hh \
  undefine-2.0.hh  undefine-3.0.hh \
  entity.hh  entity.cc  entitypointer.hh  hierarchiciterator.hh \
  algebra.hh  geometry.hh  geometry.cc  projection.hh  transformation.hh \
  leveliterator.hh  leafiterator.hh  treeiterator.hh  treeiterator.cc \
  intersection.hh  intersection.cc  intersectioniterator.hh \
  capabilities.hh  gridfactory.hh  dgfparser.hh  albertareader.hh

SUBDIRS = test

include $(top_srcdir)/am/global-rules
