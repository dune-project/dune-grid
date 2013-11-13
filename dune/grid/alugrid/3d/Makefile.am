alu3dgriddir = $(includedir)/dune/grid/alugrid/3d
alu3dgrid_HEADERS = alu3dinclude.hh alu3dgridfactory.hh alu3dgridfactory.cc \
  datahandle.hh grid.hh grid_imp.cc grid_inline.hh \
  mappings.hh mappings_imp.cc topology.hh \
  entity.hh entity_imp.cc entity_inline.hh entityseed.hh \
  faceutility.hh faceutility_imp.cc geometry.hh geometry_imp.cc \
  indexsets.hh iterator.hh iterator.cc iterator_imp.cc alu3diterators.hh \
  lbdatahandle.hh alugrid.hh capabilities.hh

headercheck_IGNORE = $(alu3dgrid_HEADERS)

if ALUGRID
  noinst_LTLIBRARIES = libalu3d.la
  libalu3d_la_SOURCES = mappings.cc topology.cc faceutility.cc

  libalu3d_la_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS)
  libalu3d_la_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS)
  libalu3d_la_LIBADD = $(ALUGRID_LIBS) $(DUNE_LIBS)
endif

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
