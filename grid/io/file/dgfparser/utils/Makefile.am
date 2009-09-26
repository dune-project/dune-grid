AM_CPPFLAGS =  $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) $(MPI_LIBS) $(LOCAL_LIBS)

noinst_PROGRAMS = dgf2dgf gmsh2dgf

dgf2dgf_SOURCES = dgf2dgf.cc
gmsh2dgf_SOURCES = gmsh2dgf.cc

include $(top_srcdir)/am/global-rules
