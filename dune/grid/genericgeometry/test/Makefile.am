GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
# GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
# GRIDTYPE = UGGRID
GRIDDIM = 3
#WORLDDIM = 3
DIMENSION = 3

check_PROGRAMS = testsubtopologies testmapping testgeo grefelement

AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS)
LDADD = $(ALL_PKG_LIBS) @LDADD@

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(AM_CPPFLAGS)			\
	-DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc

testgeo_SOURCES = testgeo.cc

grefelement_SOURCES = grefelement.cc

include $(top_srcdir)/am/global-rules
