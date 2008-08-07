GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
GRIDTYPE = YASPGRID
# GRIDTYPE = UGGRID
GRIDDIM = 3
#WORLDDIM = 3
# DIMENSION = 3

check_PROGRAMS = testsubtopologies testmapping testgeo

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc
testmapping_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
testmapping_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

testgeo_SOURCES = testgeo.cc
testgeo_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(PSGINC)
testgeo_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
