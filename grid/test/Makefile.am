# $Id$

#
## Find out which tests can be built
#

if ALBERTA
  APROG = test-alberta
endif

if ALUGRID
  ALUPROG = test-alu3dgrid test-alu3dsimplex test-alu3dcube test-alu2dsimplex
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

ALBERTA_DIM = $(GRIDDIM)

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
NORMALTESTS = test-sgrid test-oned test-yaspgrid $(APROG) $(UPROG) $(ALUPROG) $(DGFALU_UGGRID)

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

#
## common flags
#

GRIDDIM=2

# paranoia
DUNE_EXTRA_CHECKS = -DDUNE_DEVEL_MODE
# output coverage
#COVERAGE = -fprofile-arcs -ftest-coverage
AM_CXXFLAGS = $(COVERAGE)
AM_CPPFLAGS = @AM_CPPFLAGS@ $(DUNE_EXTRA_CHECKS)

#
## define the programs
#

test_sgrid_SOURCES = test-sgrid.cc
test_sgrid_LDADD = $(LOCAL_LIBS)

test_oned_SOURCES = test-oned.cc
test_oned_LDADD = $(LOCAL_LIBS)

test_yaspgrid_SOURCES = test-yaspgrid.cc
test_yaspgrid_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)
test_yaspgrid_LDFLAGS = $(AM_LDFLAGS) $(MPI_LDFLAGS)
test_yaspgrid_LDADD = $(LOCAL_LIBS) $(MPI_LIBS)

# this implicitly checks the autoconf-test as well...
test_alberta_SOURCES = test-alberta.cc
test_alberta_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) $(GRAPE_CPPFLAGS)
test_alberta_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(GRAPE_LDFLAGS)
test_alberta_LDADD = $(LOCAL_LIBS) $(ALBERTA_LIBS) $(GRAPE_LIBS)

# files for alu3dgrid
test_alu3dgrid_SOURCES = test-alu3dgrid.cc
test_alu3dgrid_CXXFLAGS = $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
test_alu3dgrid_LDADD = $(LOCAL_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(MPI_LDFLAGS)

test_alu3dsimplex_SOURCES = test-alu3dsimplex.cc
if ALBERTA
test_alu3dsimplex_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS) $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
test_alu3dsimplex_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS) $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS)
test_alu3dsimplex_LDADD = $(LOCAL_LIBS) $(ALUGRID_LIBS) $(MPI_LIBS) -lalberta_3d $(ALBERTA_BASE_LIBS)
else
test_alu3dsimplex_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS)
test_alu3dsimplex_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS) $(MPI_LDFLAGS)
test_alu3dsimplex_LDADD = $(LOCAL_LIBS) $(ALUGRID_LIBS) $(MPI_LIBS)
endif

test_alu3dcube_SOURCES = test-alu3dcube.cc
test_alu3dcube_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS)
test_alu3dcube_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS) $(MPI_LDFLAGS)
test_alu3dcube_LDADD = $(LOCAL_LIBS) $(ALUGRID_LIBS) $(MPI_LIBS)

test_alu2dsimplex_SOURCES = test-alu2dsimplex.cc
test_alu2dsimplex_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS)
test_alu2dsimplex_LDFLAGS = $(AM_LDFLAGS) $(ALUGRID_LDFLAGS) $(MPI_LDFLAGS)
test_alu2dsimplex_LDADD = $(LOCAL_LIBS) $(ALUGRID_LIBS) $(MPI_LIBS)

# libdune contains both libugX2 and libugX3, always test both dimensions
test_ug_SOURCES = test-ug.cc
test_ug_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(MPI_CPPFLAGS)
test_ug_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(MPI_LDFLAGS)
test_ug_LDADD = $(LOCAL_LIBS) $(UG_LIBS) $(MPI_LIBS)

# libdune contains both libugX2 and libugX3, always test both dimensions
test_parallel_ug_SOURCES = test-parallel-ug.cc
test_parallel_ug_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(MPI_CPPFLAGS)
test_parallel_ug_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(MPI_LDFLAGS)
test_parallel_ug_LDADD = $(LOCAL_LIBS) $(UG_LIBS) $(MPI_LIBS)

test_ug_lgm_SOURCES = test-ug-lgm.cc
test_ug_lgm_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(MPI_CPPFLAGS)
test_ug_lgm_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(MPI_LDFLAGS)
test_ug_lgm_LDADD = $(LOCAL_LIBS) $(UG_LIBS) $(MPI_LIBS)

# Test whether you can combine the different implementations in one file
test_dgfalu_uggrid_combination_SOURCES = test-dgfalu-uggrid-combination.cc
test_dgfalu_uggrid_combination_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(ALUGRID_CPPFLAGS)
test_dgfalu_uggrid_combination_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(ALUGRID_LDFLAGS)
test_dgfalu_uggrid_combination_LDADD = $(LOCAL_LIBS) $(UG_LIBS) $(ALUGRID_LIBS)

## distribution tarball
SOURCES = gridcheck.cc checkindexset.cc checkgeometryinfather.cc checkintersectionit.cc checkcommunicate.cc checkiterators.cc checktwists.cc check-albertareader.cc basicunitcube.hh checkadaptation.cc checkgeometry.cc

GRIDFILES = \
  simplex-testgrid-2-2.dgf simplex-testgrid-3-3-large.dgf simplex-testgrid-3-3.dgf \
  cube-testgrid-3-3-large.dgf

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = gridcheck.cc $(GRIDFILES) $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules
