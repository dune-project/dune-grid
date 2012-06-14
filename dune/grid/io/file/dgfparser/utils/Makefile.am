AM_CPPFLAGS = @AM_CPPFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@
LDADD = @LDADD@

EXTRA_PROGRAMS = dgf2dgf gmsh2dgf

dgf2dgf_SOURCES = dgf2dgf.cc
dgf2dgf_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DYASPGRID -DGRIDDIM=2

gmsh2dgf_SOURCES = gmsh2dgf.cc
gmsh2dgf_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DONEDGRID -DGRIDDIM=1

if ALUGRID
EXTRA_PROGRAMS += gmsh2dgfalusimplex
gmsh2dgfalusimplex_SOURCES = gmsh2dgf.cc
gmsh2dgfalusimplex_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALUGRID_CPPFLAGS)			\
	-DALUGRID_SIMPLEX -DGRIDDIM=2
gmsh2dgfalusimplex_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALUGRID_LDFLAGS)
gmsh2dgfalusimplex_LDADD =				\
	$(ALUGRID_LIBS)				\
	$(LDADD)
endif

include $(top_srcdir)/am/global-rules
