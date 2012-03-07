noinst_HEADERS =				\
	main.hh

headercheck_IGNORE =				\
	main.hh

bin_PROGRAMS =

if ALUGRID
bin_PROGRAMS += dune-structured-to-alusimplex-3d
dune_structured_to_alusimplex_3d_SOURCES = alusimplex-3d.cc
dune_structured_to_alusimplex_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_structured_to_alusimplex_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_structured_to_alusimplex_3d_LDADD =	\
	$(ALUGRID_LIBS)				\
	$(LDADD)
endif ALUGRID

include $(top_srcdir)/am/global-rules
