bin_PROGRAMS =
if EXTRA_UTILS
bin_PROGRAMS += $(EXTRAUTILS)
endif EXTRA_UTILS

EXTRA_PROGRAMS =

EXTRAUTILS =

if ALBERTA_2D
EXTRAUTILS += dune-gridinfo-gmsh-alberta-2d
endif ALBERTA_2D
EXTRA_PROGRAMS += dune-gridinfo-gmsh-alberta-2d
dune_gridinfo_gmsh_alberta_2d_SOURCES = gridinfo-gmsh-alberta-2d.cc
dune_gridinfo_gmsh_alberta_2d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA2D_CPPFLAGS)
dune_gridinfo_gmsh_alberta_2d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA2D_LDFLAGS)
dune_gridinfo_gmsh_alberta_2d_LDADD =		\
	$(ALBERTA2D_LIBS)			\
	$(LDADD)

if ALBERTA_3D
EXTRAUTILS += dune-gridinfo-gmsh-alberta-3d
endif ALBERTA_3D
EXTRA_PROGRAMS += dune-gridinfo-gmsh-alberta-3d
dune_gridinfo_gmsh_alberta_3d_SOURCES = gridinfo-gmsh-alberta-3d.cc
dune_gridinfo_gmsh_alberta_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALBERTA3D_CPPFLAGS)
dune_gridinfo_gmsh_alberta_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA3D_LDFLAGS)
dune_gridinfo_gmsh_alberta_3d_LDADD =		\
	$(ALBERTA3D_LIBS)			\
	$(LDADD)

if ALUGRID
EXTRAUTILS += dune-gridinfo-gmsh-alucube-3d
endif ALUGRID
EXTRA_PROGRAMS += dune-gridinfo-gmsh-alucube-3d
dune_gridinfo_gmsh_alucube_3d_SOURCES = gridinfo-gmsh-alucube-3d.cc
dune_gridinfo_gmsh_alucube_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_gridinfo_gmsh_alucube_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_gridinfo_gmsh_alucube_3d_LDADD =		\
	$(ALUGRID_LIBS)				\
	$(LDADD)

if ALUGRID
EXTRAUTILS += dune-gridinfo-gmsh-alusimplex-2d
endif ALUGRID
EXTRA_PROGRAMS += dune-gridinfo-gmsh-alusimplex-2d
dune_gridinfo_gmsh_alusimplex_2d_SOURCES = gridinfo-gmsh-alusimplex-2d.cc
dune_gridinfo_gmsh_alusimplex_2d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_gridinfo_gmsh_alusimplex_2d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_gridinfo_gmsh_alusimplex_2d_LDADD =	\
	$(ALUGRID_LIBS)				\
	$(LDADD)

if ALUGRID
EXTRAUTILS += dune-gridinfo-gmsh-alusimplex-3d
endif ALUGRID
EXTRA_PROGRAMS += dune-gridinfo-gmsh-alusimplex-3d
dune_gridinfo_gmsh_alusimplex_3d_SOURCES = gridinfo-gmsh-alusimplex-3d.cc
dune_gridinfo_gmsh_alusimplex_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(ALUGRID_CPPFLAGS)
dune_gridinfo_gmsh_alusimplex_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)
dune_gridinfo_gmsh_alusimplex_3d_LDADD =	\
	$(ALUGRID_LIBS)				\
	$(LDADD)

if UG
EXTRAUTILS += dune-gridinfo-gmsh-ug-2d
endif UG
EXTRA_PROGRAMS += dune-gridinfo-gmsh-ug-2d
dune_gridinfo_gmsh_ug_2d_SOURCES = gridinfo-gmsh-ug-2d.cc
dune_gridinfo_gmsh_ug_2d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(UG_CPPFLAGS)
dune_gridinfo_gmsh_ug_2d_LDFLAGS = $(AM_LDFLAGS)	\
	$(UG_LDFLAGS)
dune_gridinfo_gmsh_ug_2d_LDADD =		\
	$(UG_LIBS)				\
	$(LDADD)

if UG
EXTRAUTILS += dune-gridinfo-gmsh-ug-3d
endif UG
EXTRA_PROGRAMS += dune-gridinfo-gmsh-ug-3d
dune_gridinfo_gmsh_ug_3d_SOURCES = gridinfo-gmsh-ug-3d.cc
dune_gridinfo_gmsh_ug_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(UG_CPPFLAGS)
dune_gridinfo_gmsh_ug_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(UG_LDFLAGS)
dune_gridinfo_gmsh_ug_3d_LDADD =		\
	$(UG_LIBS)				\
	$(LDADD)

include $(top_srcdir)/am/global-rules
