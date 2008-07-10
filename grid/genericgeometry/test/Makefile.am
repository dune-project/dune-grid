GEOMETRYTYPE = cube
# GRIDTYPE = ALBERTAGRID
GRIDTYPE = YASPGRID
# GRIDTYPE = PARALLELSIMPLEXGRID
GRIDDIM = 3
DIMENSION = 3

PSGINC = -I../../../../dune-psg/

check_PROGRAMS = testsubtopologies testmapping

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmapping_SOURCES = testmapping.cc
testmapping_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(PSGINC)
testmapping_LDADD = $(LOCAL_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)

