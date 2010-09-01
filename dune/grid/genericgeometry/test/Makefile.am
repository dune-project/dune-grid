GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
# GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
# GRIDTYPE = UGGRID
GRIDDIM = 3
#WORLDDIM = 3
DIMENSION = 3

check_HEADERS = testgeo.hh

check_PROGRAMS = testsubtopologies testmapping testgeo testmatrix
TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS)
LDADD = $(ALL_PKG_LIBS) @LDADD@

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(AM_CPPFLAGS)			\
	-DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc testgeo.hh

testgeo_SOURCES = testgeo.cc testgeo.hh

testmatrix_SOURCES = testmatrix.cc
testmatrix_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
testmatrix_LDADD = $(LDADD) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules
