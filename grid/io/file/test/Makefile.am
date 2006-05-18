# $Id$

#
## estimate which tests can be built
#

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
NORMALTESTS =

# list of tests to run (yaspgrid is special case)
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

#
## common flags
#

# paranoia
DUNE_EXTRA_CHECKS =  -DDUNE_ISTL_WITH_CHECKING -DDUNE_DEVEL_MODE

AM_CPPFLAGS = @AM_CPPFLAGS@ $(DUNE_EXTRA_CHECKS)

#
## define the programs
#

#
## distribution tarball

testdir = $(includedir)/dune/grid/io/file/test/
test_HEADERS =

GRIDFILES = ug-testgrid-2.am ug-testgrid-3.am

# gridcheck not used explicitly, we should still ship it :)
#EXTRA_DIST = gridcheck.cc $(GRIDFILES)

CLEANFILES = *.gcda *.gcno semantic.cache

include $(top_srcdir)/am/global-rules
