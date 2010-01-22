# $Id$

if UG
TESTPROGS = mcmgmappertest
endif

# which tests to run
TESTS = test-quadrature test-genericquadrature referenceelementstest checkrefnumbering scsgmappertest $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CXXFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

mcmgmappertest_SOURCES = mcmgmappertest.cc
mcmgmappertest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(UG_CPPFLAGS)				\
	$(AMIRAMESH_CPPFLAGS)
mcmgmappertest_LDFLAGS = $(AM_LDFLAGS)		\
	$(UG_LDFLAGS)				\
	$(AMIRAMESH_LDFLAGS)
mcmgmappertest_LDADD =				\
	$(AMIRAMESH_LIBS)			\
	$(UG_LIBS)				\
	$(LDADD)

scsgmappertest_SOURCES = scsgmappertest.cc

referenceelementstest_SOURCES = referenceelementstest.cc

checkrefnumbering_SOURCES = checkrefnumbering.cc

test_quadrature_SOURCES = test-quadrature.cc

test_genericquadrature_SOURCES = test-genericquadrature.cc

include $(top_srcdir)/am/global-rules
