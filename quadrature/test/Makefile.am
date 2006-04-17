# $Id$

#
## define the lists of tests to build and run
#

TESTS = test-quadrature

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

#
## common flags
#

#
## define the programs
#

test_quadrature_SOURCES = test-quadrature.cc
test_quadrature_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
