# $Id$

if UG
TESTPROGS = mcmgmappertest
endif

# which tests to run
TESTS = test-quadrature test-genericquadrature referenceelementstest checkrefnumbering scsgmappertest $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

mcmgmappertest_SOURCES = mcmgmappertest.cc
mcmgmappertest_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
mcmgmappertest_LDADD = $(LOCAL_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS)

scsgmappertest_SOURCES = scsgmappertest.cc
scsgmappertest_LDADD = $(LOCAL_LIBS)

referenceelementstest_SOURCES = referenceelementstest.cc
referenceelementstest_LDADD = $(LOCAL_LIBS)

checkrefnumbering_SOURCES = checkrefnumbering.cc
checkrefnumbering_LDADD = $(LOCAL_LIBS)

test_quadrature_SOURCES = test-quadrature.cc
test_quadrature_LDADD = $(LOCAL_LIBS)

test_genericquadrature_SOURCES = test-genericquadrature.cc
test_genericquadrature_LDADD = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
