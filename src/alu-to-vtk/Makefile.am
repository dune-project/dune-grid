noinst_HEADERS =				\
	main-volume.hh

headercheck_IGNORE =				\
	main-volume.hh

bin_PROGRAMS =

if ALUGRID
bin_PROGRAMS += dune-alusimplex-3d-to-vtk-volume
bin_PROGRAMS += dune-alucube-3d-to-vtk-volume
dune_alusimplex_3d_to_vtk_volume_SOURCES = alusimplex-3d-volume.cc
dune_alusimplex_3d_to_vtk_volume_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_alusimplex_3d_to_vtk_volume_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_alusimplex_3d_to_vtk_volume_LDADD =	\
	$(ALUGRID_LIBS)				\
	$(LDADD)
dune_alucube_3d_to_vtk_volume_SOURCES = alucube-3d-volume.cc
dune_alucube_3d_to_vtk_volume_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_alucube_3d_to_vtk_volume_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_alucube_3d_to_vtk_volume_LDADD =	\
	$(ALUGRID_LIBS)				\
	$(LDADD)
endif ALUGRID

include $(top_srcdir)/am/global-rules
