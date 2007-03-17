# $Id$

#
## estimate which tests can be built
#

if MPI
  YPROG = test-yaspgrid
endif

if ALBERTA
  APROG = test-alberta
endif

if ALUGRID
  ALUPROG = test-alu3dgrid
endif

if UG
if UG_LGMDOMAIN
  UPROG = test-ug-lgm
else
  UPROG = test-ug
endif
if ALUGRID
  DGFALU_UGGRID = test-dgfalu-uggrid-combination
else
  DGFALU_UGGRID =
endif
endif

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
NORMALTESTS = test-sgrid test-oned $(APROG) $(UPROG) $(ALUPROG) $(DGFALU_UGGRID)

# list of tests to run (yaspgrid is special case)
TESTS = $(NORMALTESTS) $(YPROG)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS) $(YPROG)

#
## common flags
#

# paranoia
DUNE_EXTRA_CHECKS = -DDUNE_DEVEL_MODE
# output coverage
#COVERAGE = -fprofile-arcs -ftest-coverage
AM_CPPFLAGS = @AM_CPPFLAGS@ $(COVERAGE) $(DUNE_EXTRA_CHECKS)

#
## define the programs
#

test_sgrid_SOURCES = test-sgrid.cc
test_sgrid_DEPENDENCIES = $(LOCAL_LIBS)
test_sgrid_LDADD = $(LOCAL_LIBS)

test_oned_SOURCES = test-oned.cc
test_oned_DEPENDENCIES = $(LOCAL_LIBS)
test_oned_LDADD = $(LOCAL_LIBS)

test_yaspgrid_SOURCES = test-yaspgrid.cc
test_yaspgrid_CXXFLAGS = $(AM_CPPFLAGS) $(MPI_CPPFLAGS)
test_yaspgrid_LDADD = $(LOCAL_LIBS) $(MPI_LDFLAGS)
test_yaspgrid_DEPENDENCIES = $(LOCAL_LIBS)

# this implicitly checks the autoconf-test as well...
test_alberta_SOURCES = test-alberta.cc
test_alberta_CXXFLAGS = $(ALBERTA_CPPFLAGS) $(AM_CPPFLAGS)
test_alberta_LDADD = $(LOCAL_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
test_alberta_DEPENDENCIES = $(LOCAL_LIBS)

# files for alu3dgrid
test_alu3dgrid_SOURCES = test-alu3dgrid.cc
test_alu3dgrid_CXXFLAGS = $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
test_alu3dgrid_LDADD = $(LOCAL_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(MPI_LDFLAGS)
test_alu3dgrid_DEPENDENCIES = $(LOCAL_LIBS)

# libdune contains both libugX2 and libugX3, always test both dimensions
test_ug_SOURCES = test-ug.cc
test_ug_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(MPI_CPPFLAGS)
test_ug_LDADD = $(LOCAL_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS)
test_ug_DEPENDENCIES = $(LOCAL_LIBS)

test_ug_lgm_SOURCES = test-ug-lgm.cc
test_ug_lgm_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(MPI_CPPFLAGS)
test_ug_lgm_LDADD = $(LOCAL_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS)
test_ug_lgm_DEPENDENCIES = $(LOCAL_LIBS)

# Test whether you can combine the different implementations in one file
test_dgfalu_uggrid_combination_SOURCES = test-dgfalu-uggrid-combination.cc
test_dgfalu_uggrid_combination_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(ALUGRID_CPPFLAGS)
test_dgfalu_uggrid_combination_LDADD = $(LOCAL_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)
test_dgfalu_uggrid_combination_DEPENDENCIES = $(LOCAL_LIBS)

## distribution tarball
SOURCES = gridcheck.cc checkindexset.cc checkgeometryinfather.cc checkintersectionit.cc checkcommunicate.cc

GRIDFILES = alberta-testgrid-3-3.al \
  simplex-testgrid-2-2.dgf simplex-testgrid-3-3-large.dgf simplex-testgrid-3-3.dgf \
  cube-testgrid-3-3-large.dgf cube-testgrid-3-3.dgf

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = gridcheck.cc $(GRIDFILES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log defaults

include $(top_srcdir)/am/global-rules
