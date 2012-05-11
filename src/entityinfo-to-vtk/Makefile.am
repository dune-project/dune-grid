noinst_HEADERS =				\
	entityinfo-to-vtk.hh			\
	entityinfo-to-vtk-main-yasp.hh

bin_PROGRAMS =

bin_PROGRAMS += dune-entityinfo-to-vtk-yasp-3d
dune_entityinfo_to_vtk_yasp_3d_SOURCES = dune-entityinfo-to-vtk-yasp-3d.cc
dune_entityinfo_to_vtk_yasp_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)
dune_entityinfo_to_vtk_yasp_3d_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)
dune_entityinfo_to_vtk_yasp_3d_LDADD =		\
	$(DUNEMPILIBS)				\
	$(LDADD)

include $(top_srcdir)/am/global-rules
