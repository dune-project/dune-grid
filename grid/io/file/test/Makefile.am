ALLTESTS = vtktest gnuplottest

# Currently, Star-CD files can only be read into UGGrid objects
if UG
ALLTESTS += starcdreadertest
endif

AM_CXXFLAGS= -g -Wall -O0

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

starcdreadertest_SOURCES = starcdreadertest.cc
starcdreadertest_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)
starcdreadertest_LDADD = $(LOCAL_LIBS)
starcdreadertest_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp
