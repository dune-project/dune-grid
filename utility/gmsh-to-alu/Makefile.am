noinst_HEADERS =				\
	main.hh

bin_PROGRAMS =
if EXTRA_UTILS
bin_PROGRAMS += $(EXTRAUTILS)
endif EXTRA_UTILS

EXTRA_PROGRAMS =
EXTRAUTILS =

if ALUGRID
EXTRAUTILS += dune-gmsh-to-alusimplex-3d
endif ALUGRID
EXTRA_PROGRAMS += dune-gmsh-to-alusimplex-3d
dune_gmsh_to_alusimplex_3d_SOURCES = alusimplex-3d.cc
dune_gmsh_to_alusimplex_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_gmsh_to_alusimplex_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_gmsh_to_alusimplex_3d_LDADD =		\
	$(ALUGRID_LIBS)				\
	$(LDADD)

include $(top_srcdir)/am/global-rules
