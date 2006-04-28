# $Id$

if UG
if AMIRAMESH
TESTPROGS = mcmgmappertest
endif
endif

# which tests to run
TESTS = test-quadrature $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

mcmgmappertest_SOURCES = mcmgmappertest.cc
mcmgmappertest_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
mcmgmappertest_LDADD = $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS)
mcmgmappertest_DEPENDENCIES = $(LOCAL_LIBS)

test_quadrature_SOURCES = test-quadrature.cc
test_quadrature_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
