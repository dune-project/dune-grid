GEOMETRYTYPE = cube
GRIDTYPE = ALBERTAGRID
# GRIDTYPE = YASPGRID
GRIDDIM = 3
DIMENSION = 3

check_PROGRAMS = testsubtopologies testmapping

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc
testmapping_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
testmapping_LDADD = $(LOCAL_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

