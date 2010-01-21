# $Id$

if ALBERTA
noinst_LTLIBRARIES = libalbertagrid_1d.la libalbertagrid_2d.la libalbertagrid_3d.la

libalbertagrid_sources = dgfparser.cc  elementinfo.cc  indexsets.cc  macrodata.cc  meshpointer.cc
libalbertagrid_cppflags = $(AM_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DENABLE_ALBERTA

libalbertagrid_1d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_1d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=1
libalbertagrid_1d_la_LIBADD = $(DUNE_COMMON_LIBS) -L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_1d

libalbertagrid_2d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_2d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=2
libalbertagrid_2d_la_LIBADD = $(DUNE_COMMON_LIBS) -L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_2d

libalbertagrid_3d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_3d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=3
libalbertagrid_3d_la_LIBADD = $(DUNE_COMMON_LIBS) -L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_3d
endif

albertadir = $(includedir)/dune/grid/albertagrid/
alberta_HEADERS = agrid.hh albertagrid.cc \
  agmemory.hh  albertaheader.hh \
  indexsets.hh datahandle.hh \
  misc.hh  macroelement.hh  elementinfo.hh  geometrycache.hh  meshpointer.hh \
  macrodata.hh  dofadmin.hh  dofvector.hh  refinement.hh  coordcache.hh \
  level.hh \
  undefine-2.0.hh  undefine-3.0.hh \
  entity.hh  entity.cc  entitypointer.hh  hierarchiciterator.hh \
  algebra.hh  geometry.hh  geometry.cc  projection.hh  transformation.hh \
  leveliterator.hh  leafiterator.hh  treeiterator.hh \
  intersection.hh  intersection.cc  intersectioniterator.hh \
  capabilities.hh  gridfactory.hh  dgfparser.hh  albertareader.hh \
  gridfamily.hh

SUBDIRS = test

include $(top_srcdir)/am/global-rules
