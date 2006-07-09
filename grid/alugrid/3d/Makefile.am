AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALUGRID_CPPFLAGS)

if ALUGRID
  noinst_LTLIBRARIES = libalu3d.la
  libalu3d_la_SOURCES = mappings.cc topology.cc faceutility.cc
  #libalu3d_la_LIBADD = -lcommon
else
  noinst_LTLIBRARIES = libalu3d.la
  libalu3d_la_SOURCES = topology.cc
  #libalu3d_la_LIBADD = -lcommon
endif

alu3dgriddir = $(includedir)/dune/grid/alugrid/3d
alu3dgrid_HEADERS = alu3dinclude.hh \
  datahandle.hh memory.hh grid.hh grid_imp.cc \
  mappings.hh mappings_imp.cc topology.hh entity.hh entity_imp.cc \
  faceutility.hh faceutility_imp.cc geometry.hh geometry_imp.cc \
  myautoptr.hh indexsets.hh iterator.hh iterator_imp.cc alu3diterators.hh \
  alu3dutility.hh

include $(top_srcdir)/am/global-rules
