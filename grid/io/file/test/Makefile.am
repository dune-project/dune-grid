ALLTESTS = vtktest

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS)

vtktest_SOURCES = vtktest.cc
vtktest_LDADD = $(LOCAL_LIBS)
vtktest_DEPENDENCIES = $(LOCAL_LIBS)

CLEANFILES=*.vtu *.vtp

include $(top_srcdir)/am/global-rules
