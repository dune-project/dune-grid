AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(ALL_PKG_CPPFLAGS)			\
	$(DUNEMPICPPFLAGS)
LDADD = $(LOCAL_LIBS)				\
	$(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)	\
	$(DUNEMPILIBS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)			\
	$(ALL_PKG_LDFLAGS)

#don't build as noinst_PROGRAMS because they depend on the lib
#noinst_PROGRAMS = dgf2dgf gmsh2dgf
check_PROGRAMS = dgf2dgf gmsh2dgf

dgf2dgf_SOURCES = dgf2dgf.cc
gmsh2dgf_SOURCES = gmsh2dgf.cc

include $(top_srcdir)/am/global-rules
