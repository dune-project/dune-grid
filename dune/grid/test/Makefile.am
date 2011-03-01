# $Id$

#
## Find out which tests can be built
#

if ALBERTA
  APROG = test-alberta test-alberta-generic
endif

if ALUGRID
  ALUPROG = test-alugrid
endif

if UG
if UG_LGMDOMAIN
  UPROG = test-ug-lgm
else
  UPROG = test-ug test-parallel-ug
endif
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
NORMALTESTS = test-sgrid test-oned test-yaspgrid test-geogrid $(APROG) $(UPROG) $(ALUPROG) $(DGFALU_UGGRID)

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

#
## common flags
#

# paranoia
DUNE_EXTRA_CHECKS = -DDUNE_DEVEL_MODE
# output coverage
#COVERAGE = -fprofile-arcs -ftest-coverage
AM_CXXFLAGS = $(COVERAGE)
AM_CPPFLAGS = @AM_CPPFLAGS@ $(DUNE_EXTRA_CHECKS)

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

#
## define the programs
#

test_sgrid_SOURCES = test-sgrid.cc

test_oned_SOURCES = test-oned.cc

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
test_alberta_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) $(GRAPE_CPPFLAGS)
test_alberta_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(GRAPE_LDFLAGS)
test_alberta_LDADD =				\
	$(GRAPE_LIBS)				\
	$(ALBERTA_LIBS)				\
	$(LDADD)

test_alberta_generic_SOURCES = $(test_alberta_SOURCES)
test_alberta_generic_CPPFLAGS = $(test_alberta_CPPFLAGS) -DDUNE_ALBERTA_USE_GENERICGEOMETRY=1
test_alberta_generic_LDFLAGS = $(test_alberta_LDFLAGS)
test_alberta_generic_LDADD = $(test_alberta_LDADD)

# files for alugrid
test_alugrid_SOURCES = test-alugrid.cc
test_alugrid_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALL_PKG_CPPFLAGS)
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
	-DCACHECOORDFUNCTION=$(CACHECOORDFUNCTION)
test_geogrid_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALL_PKG_LDFLAGS)
test_geogrid_LDADD =				\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

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

test_ug_lgm_SOURCES = test-ug-lgm.cc
test_ug_lgm_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)			\
	$(UG_CPPFLAGS)
test_ug_lgm_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(UG_LDFLAGS)
test_ug_lgm_LDADD =				\
	$(UG_LIBS)				\
	$(DUNEMPILIBS)				\
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
SOURCES = gridcheck.cc staticcheck.hh checkindexset.cc checkgeometryinfather.cc checkintersectionit.cc checkcommunicate.cc checkiterators.cc checktwists.cc check-albertareader.cc basicunitcube.hh checkadaptation.cc checkgeometry.cc checkpartition.cc

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules
