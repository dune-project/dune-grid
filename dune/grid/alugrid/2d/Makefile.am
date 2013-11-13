alu2dgriddir = $(includedir)/dune/grid/alugrid/2d
alu2dgrid_HEADERS = alu2dinclude.hh  entity_imp.cc  geometry_imp.cc \
  grid_imp.cc  intersection.hh  iterator.hh  entity.hh  geometry.hh   grid.hh \
  indexsets.hh  intersection_imp.cc  iterator_imp.cc datahandle.hh \
  alu2dgridfactory.hh  alugrid.hh capabilities.hh  bndprojection.hh \
  intersection.cc  alu2dgridfactory.cc  entityseed.hh

headercheck_IGNORE = $(alu2dgrid_HEADERS)

if ALUGRID
  noinst_LTLIBRARIES = libalu2d.la
  libalu2d_la_SOURCES = intersection.cc  alu2dgridfactory.cc

  libalu2d_la_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS)
  libalu2d_la_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS)
  libalu2d_la_LIBADD = $(ALUGRID_LIBS) $(DUNE_LIBS)
endif

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
