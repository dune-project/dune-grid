GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
# GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
# GRIDTYPE = UGGRID
GRIDDIM = 3
#WORLDDIM = 3
DIMENSION = 3

check_HEADERS = testgeo.hh

check_PROGRAMS = testmapping testgeo
TESTS = testmapping testgeo

AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS) $(DUNEMPILIBS)
LDADD = $(ALL_PKG_LIBS) @LDADD@ $(DUNEMPILDFLAGS)

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

testmapping_SOURCES = testmapping.cc testgeo.hh

testgeo_SOURCES = testgeo.cc testgeo.hh

include $(top_srcdir)/am/global-rules
