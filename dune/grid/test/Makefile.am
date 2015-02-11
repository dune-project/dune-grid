#
## Find out which tests can be built
#

if ALBERTA
  APROG = test-alberta-1-1 test-alberta-1-2 test-alberta-2-2 \
          test-alberta-1-3 test-alberta-2-3 test-alberta-3-3 \
          test-alberta-generic
  ALBERTA_EXTRA_PROGS = test-alberta
endif

if ALUGRID
  ALUPROG = test-alugrid
endif

if UG
  UPROG = test-ug test-parallel-ug test-loadbalancing
if ALUGRID
  DGFALU_UGGRID = test-dgfalu-uggrid-combination
else
  DGFALU_UGGRID =
endif
endif

#
## Defines for gridtype.hh
#
GRIDTYPE=YASPGRID
GRIDDIM=2

ALBERTA_DIM = $(WORLDDIM)

#
## Some parameters for the geogrid test
#
COORDFUNCTION = IdenticalCoordFunction\<double,$(GRIDDIM)\>
# COORDFUNCTION = Helix
# COORDFUNCTION = ThickHelix
# COORDFUNCTION = DGFCoordFunction\<$(GRIDDIM),$(GRIDDIM)\>
# COORDFUNCTION = DGFCoordFunction\<2,3\>
CACHECOORDFUNCTION = 0

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
NORMALTESTS = test-sgrid test-oned test-yaspgrid test-geogrid \
              test-identitygrid test-mcmg-geogrid \
              $(APROG) $(UPROG) $(ALUPROG) $(DGFALU_UGGRID)

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

EXTRA_PROGRAMS = $(ALBERTA_EXTRA_PROGS)

#
## common flags
#

# output coverage
#COVERAGE = -fprofile-arcs -ftest-coverage
AM_CXXFLAGS = $(COVERAGE)
AM_CPPFLAGS = @AM_CPPFLAGS@ -DHAVE_DUNE_GRID=1

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

#
## define the programs
#

test_sgrid_SOURCES = test-sgrid.cc

test_oned_SOURCES = test-oned.cc

test_identitygrid_SOURCES = test-identitygrid.cc

test_yaspgrid_SOURCES = test-yaspgrid.cc
test_yaspgrid_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)
test_yaspgrid_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)
test_yaspgrid_LDADD =				\
	$(DUNEMPILIBS)				\
	$(LDADD)

# this implicitly checks the autoconf-test as well...
test_alberta_SOURCES = test-alberta.cc
test_alberta_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) $(GRAPE_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(GRAPE_LDFLAGS)
test_alberta_LDADD =				\
	$(GRAPE_LIBS)				\
	$(ALBERTA_LIBS)				\
	$(LDADD)

test_alberta_1_1_SOURCES = test-alberta.cc
test_alberta_1_1_CPPFLAGS = -DGRIDDIM=1 -DWORLDDIM=1 $(ALBERTA1D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_1_1_LDFLAGS = $(ALBERTA1D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_1_1_LDADD = $(ALBERTA1D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_1_2_SOURCES = test-alberta.cc
test_alberta_1_2_CPPFLAGS = -DGRIDDIM=1 -DWORLDDIM=2 $(ALBERTA2D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_1_2_LDFLAGS = $(ALBERTA2D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_1_2_LDADD = $(ALBERTA2D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_2_2_SOURCES = test-alberta.cc
test_alberta_2_2_CPPFLAGS = -DGRIDDIM=2 -DWORLDDIM=2 $(ALBERTA2D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_2_2_LDFLAGS = $(ALBERTA2D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_2_2_LDADD = $(ALBERTA2D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_1_3_SOURCES = test-alberta.cc
test_alberta_1_3_CPPFLAGS = -DGRIDDIM=1 -DWORLDDIM=3 $(ALBERTA3D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_1_3_LDFLAGS = $(ALBERTA3D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_1_3_LDADD = $(ALBERTA3D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_2_3_SOURCES = test-alberta.cc
test_alberta_2_3_CPPFLAGS = -DGRIDDIM=2 -DWORLDDIM=3 $(ALBERTA3D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_2_3_LDFLAGS = $(ALBERTA3D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_2_3_LDADD = $(ALBERTA3D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_3_3_SOURCES = test-alberta.cc
test_alberta_3_3_CPPFLAGS = -DGRIDDIM=3 -DWORLDDIM=3 $(ALBERTA3D_CPPFLAGS) $(GRAPE_CPPFLAGS) $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_3_3_LDFLAGS = $(ALBERTA3D_LDFLAGS) $(GRAPE_LDFLAGS) $(AM_LDFLAGS)
test_alberta_3_3_LDADD = $(ALBERTA3D_LIBS) $(GRAPE_LIBS) $(LDADD)

test_alberta_generic_SOURCES = $(test_alberta_SOURCES)
test_alberta_generic_CPPFLAGS = $(test_alberta_CPPFLAGS) -DDUNE_ALBERTA_USE_GENERICGEOMETRY=1 \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alberta_generic_LDFLAGS = $(test_alberta_LDFLAGS)
test_alberta_generic_LDADD = $(test_alberta_LDADD)

# files for alugrid
test_alugrid_SOURCES = test-alugrid.cc
test_alugrid_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALL_PKG_CPPFLAGS)                     \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_alugrid_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALL_PKG_LDFLAGS)			\
	-static
test_alugrid_LDADD =				\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

test_geogrid_SOURCES = test-geogrid.cc functions.hh
test_geogrid_CPPFLAGS = $(AM_CPPFLAGS)			\
	$(ALL_PKG_CPPFLAGS)				\
	-DCOORDFUNCTION=$(COORDFUNCTION)		\
	-DCACHECOORDFUNCTION=$(CACHECOORDFUNCTION)      \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_geogrid_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALL_PKG_LDFLAGS)
test_geogrid_LDADD =				\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

test_mcmg_geogrid_SOURCES = test-mcmg-geogrid.cc
test_mcmg_geogrid_CPPFLAGS = $(AM_CPPFLAGS) \
	-DDUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
test_mcmg_geogrid_LDFLAGS = $(AM_LDFLAGS)
test_mcmg_geogrid_LDADD = $(LDADD)

# libdune contains both libugX2 and libugX3, always test both dimensions
test_ug_SOURCES = test-ug.cc
test_ug_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)			\
	$(UG_CPPFLAGS)
test_ug_LDFLAGS = $(AM_LDFLAGS)			\
	$(DUNEMPILDFLAGS)			\
	$(UG_LDFLAGS)
test_ug_LDADD =					\
	$(UG_LIBS)				\
	$(DUNEMPILIBS)				\
	$(LDADD)

# libdune contains both libugX2 and libugX3, always test both dimensions
test_parallel_ug_SOURCES = test-parallel-ug.cc
test_parallel_ug_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(UG_CPPFLAGS)
test_parallel_ug_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)			\
	$(UG_LDFLAGS)
test_parallel_ug_LDADD =			\
	$(UG_LIBS)				\
	$(DUNEMPILIBS)				\
	$(LDADD)

test_loadbalancing_SOURCES = test-loadbalancing.cc
test_loadbalancing_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(UG_CPPFLAGS)				\
	$(PARMETIS_CPPFLAGS)
test_loadbalancing_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)			\
	$(UG_LDFLAGS)				\
	$(PARMETIS_LDFLAGS)
test_loadbalancing_LDADD =			\
	$(UG_LIBS)				\
	$(DUNEMPILIBS)				\
	$(PARMETIS_LIBS)			\
	$(LDADD)

# Test whether you can combine the different implementations in one file
test_dgfalu_uggrid_combination_SOURCES = test-dgfalu-uggrid-combination.cc
test_dgfalu_uggrid_combination_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(ALUGRID_CPPFLAGS)
test_dgfalu_uggrid_combination_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(ALUGRID_LDFLAGS)
test_dgfalu_uggrid_combination_LDADD =		\
	$(ALUGRID_LIBS)				\
	$(UG_LIBS)				\
	$(LDADD)

## distribution tarball
SOURCES = basicunitcube.hh                      \
          check-albertareader.cc                \
          check-albertareader.hh                \
          checkadaptation.cc                    \
          checkadaptation.hh                    \
          checkcommunicate.cc                   \
          checkcommunicate.hh                   \
          checkentitylifetime.hh                \
          checkentityseed.cc                    \
          checkentityseed.hh                    \
          checkgeometry.cc                      \
          checkgeometry.hh                      \
          checkgeometryinfather.cc              \
          checkgeometryinfather.hh              \
          checkindexset.cc                      \
          checkindexset.hh                      \
          checkintersectionit.cc                \
          checkintersectionit.hh                \
          checkintersectionlifetime.hh          \
          checkiterators.cc                     \
          checkiterators.hh                     \
          checkjacobians.cc                     \
          checkjacobians.hh                     \
          checkpartition.cc                     \
          checkpartition.hh                     \
          checktwists.cc                        \
          checktwists.hh                        \
          functions.hh                          \
          gridcheck.cc                          \
          gridcheck.hh                          \
          staticcheck.hh

NORMALTESTS += testiteratorranges
testiteratorranges_SOURCES = testiteratorranges.cc
testiteratorranges_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)
testiteratorranges_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)
testiteratorranges_LDADD =			\
	$(DUNEMPILIBS)				\
	$(LDADD)

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = CMakeLists.txt $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

# install the tests as we want to support testing 3rdparty grids with installed dune-grid
testgriddir = $(includedir)/dune/grid/test
testgrid_HEADERS = $(SOURCES)

include $(top_srcdir)/am/global-rules
