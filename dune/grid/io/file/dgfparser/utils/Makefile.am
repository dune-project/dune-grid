AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(ALL_PKG_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(ALL_PKG_LDFLAGS)
LDADD =						\
	$(ALL_PKG_LIBS)				\
	@LDADD@

GRIDDIM = 2
GRIDTYPE = YASPGRID

#don't build as noinst_PROGRAMS because they depend on the lib
#noinst_PROGRAMS = dgf2dgf gmsh2dgf
check_PROGRAMS = dgf2dgf gmsh2dgf

dgf2dgf_SOURCES = dgf2dgf.cc
gmsh2dgf_SOURCES = gmsh2dgf.cc

include $(top_srcdir)/am/global-rules
