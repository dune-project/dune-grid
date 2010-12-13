# $Id$

if ALBERTA
libalbertagrid_sources = dgfparser.cc  elementinfo.cc  indexsets.cc  macrodata.cc  meshpointer.cc
libalbertagrid_cppflags = $(AM_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DENABLE_ALBERTA
endif

ALBERTAGRID_LIBS =

if ALBERTA_1D
ALBERTAGRID_LIBS += libalbertagrid_1d.la
endif
if ALBERTA_2D
ALBERTAGRID_LIBS += libalbertagrid_2d.la
endif
if ALBERTA_3D
ALBERTAGRID_LIBS += libalbertagrid_3d.la
endif
if ALBERTA_4D
ALBERTAGRID_LIBS += libalbertagrid_4d.la
endif
if ALBERTA_5D
ALBERTAGRID_LIBS += libalbertagrid_5d.la
endif
if ALBERTA_6D
ALBERTAGRID_LIBS += libalbertagrid_6d.la
endif
if ALBERTA_7D
ALBERTAGRID_LIBS += libalbertagrid_7d.la
endif
if ALBERTA_8D
ALBERTAGRID_LIBS += libalbertagrid_8d.la
endif
if ALBERTA_9D
ALBERTAGRID_LIBS += libalbertagrid_9d.la
endif

noinst_LTLIBRARIES = $(ALBERTAGRID_LIBS)

albertadir = $(includedir)/dune/grid/albertagrid/
alberta_HEADERS = agrid.hh albertagrid.cc \
  albertaheader.hh \
  indexsets.hh datahandle.hh \
  misc.hh  macroelement.hh  elementinfo.hh  geometrycache.hh  meshpointer.hh \
  macrodata.hh  dofadmin.hh  dofvector.hh  refinement.hh  coordcache.hh \
  level.hh \
  undefine-2.0.hh  undefine-3.0.hh \
  entity.hh  entity.cc  entitypointer.hh  entityseed.hh  hierarchiciterator.hh \
  algebra.hh  geometry.hh  geometry.cc  projection.hh  transformation.hh \
  leveliterator.hh  leafiterator.hh  treeiterator.hh \
  intersection.hh  intersection.cc  intersectioniterator.hh \
  capabilities.hh  gridfactory.hh  dgfparser.hh  albertareader.hh \
  gridfamily.hh  persistentcontainer.hh

SUBDIRS = test

libalbertagrid_1d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_1d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=1
libalbertagrid_1d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_1d $(DUNE_LIBS)

libalbertagrid_2d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_2d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=2
libalbertagrid_2d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_2d $(DUNE_LIBS)

libalbertagrid_3d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_3d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=3
libalbertagrid_3d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_3d $(DUNE_LIBS)

libalbertagrid_4d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_4d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=4
libalbertagrid_4d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_4d $(DUNE_LIBS)

libalbertagrid_5d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_5d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=5
libalbertagrid_5d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_5d $(DUNE_LIBS)

libalbertagrid_6d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_6d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=6
libalbertagrid_6d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_6d $(DUNE_LIBS)

libalbertagrid_7d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_7d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=7
libalbertagrid_7d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_7d $(DUNE_LIBS)

libalbertagrid_8d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_8d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=8
libalbertagrid_8d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_8d $(DUNE_LIBS)

libalbertagrid_9d_la_SOURCES = $(libalbertagrid_sources)
libalbertagrid_9d_la_CPPFLAGS = $(libalbertagrid_cppflags) -DALBERTA_DIM=9
libalbertagrid_9d_la_LIBADD =-L$(ALBERTAROOT)/lib $(ALBERTA_BASE_LIBS) -lalberta_9d $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules
