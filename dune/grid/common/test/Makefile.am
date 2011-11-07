# $Id$

if UG
TESTPROGS = mcmgmappertest
endif

# which tests to run
TESTS = scsgmappertest $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CXXFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

mcmgmappertest_SOURCES = mcmgmappertest.cc
mcmgmappertest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(UG_CPPFLAGS)				\
	$(AMIRAMESH_CPPFLAGS) $(DUNEMPICPPFLAGS)
mcmgmappertest_LDFLAGS = $(AM_LDFLAGS)		\
	$(UG_LDFLAGS)				\
	$(AMIRAMESH_LDFLAGS)
mcmgmappertest_LDADD =				\
	$(AMIRAMESH_LIBS)			\
	$(UG_LIBS)				\
	$(LDADD)

scsgmappertest_SOURCES = scsgmappertest.cc

include $(top_srcdir)/am/global-rules
