libalu3d_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
libalu3d_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALUGRID_LDFLAGS)
libalu3d_la_LIBADD = 				\
	$(ALUGRID_LIBS)				\
	$(DUNE_LIBS)

noinst_LTLIBRARIES = libalu3d.la
if ALUGRID
  libalu3d_la_SOURCES = mappings.cc mappings_imp.cc topology.cc faceutility.cc \
                        alu3dgridfactory.cc grid_imp.cc
else
  libalu3d_la_SOURCES = topology.cc
endif

alu3dgriddir = $(includedir)/dune/grid/alugrid/3d
alu3dgrid_HEADERS = alu3dinclude.hh alu3dgridfactory.hh alu3dgridfactory.cc \
  datahandle.hh memory.hh grid.hh grid_imp.cc grid_inline.hh \
  mappings.hh mappings_imp.cc topology.hh entity.hh entity_imp.cc \
  faceutility.hh faceutility_imp.cc geometry.hh geometry_imp.cc \
  indexsets.hh iterator.hh iterator_imp.cc alu3diterators.hh \
  lbdatahandle.hh alugrid.hh capabilities.hh

headercheck_IGNORE = $(alu3dgrid_HEADERS)

include $(top_srcdir)/am/global-rules
