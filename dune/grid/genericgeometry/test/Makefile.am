GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
# GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
# GRIDTYPE = UGGRID
GRIDDIM = 3
#WORLDDIM = 3
DIMENSION = 3

check_HEADERS = testgeo.hh

check_PROGRAMS = testbasicgeometry testsubtopologies testmapping testgeo testmatrix
TESTS = testbasicgeometry testsubtopologies testmapping testgeo

AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS) $(DUNEMPILIBS)
LDADD = $(ALL_PKG_LIBS) @LDADD@ $(DUNEMPILDFLAGS)

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

testbasicgeometry_SOURCES = testbasicgeometry.cc

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(AM_CPPFLAGS)			\
	-DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc testgeo.hh

testgeo_SOURCES = testgeo.cc testgeo.hh

testmatrix_SOURCES = testmatrix.cc
testmatrix_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
testmatrix_LDADD = $(LDADD) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules
