ALLTESTS = vtktest gnuplottest

CXXFLAGS= -g -Wall -O0

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS)

vtktest_SOURCES = vtktest.cc
vtktest_LDADD = $(LOCAL_LIBS)
vtktest_DEPENDENCIES = $(LOCAL_LIBS)

gnuplottest_SOURCES = gnuplottest.cc
gnuplottest_LDADD = $(LOCAL_LIBS)
gnuplottest_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp
