ALLTESTS = vtktest gnuplottest vtksequencetest subsamplingvtktest gmshtest gmshtest-all

# run amiramesh test only if amiramesh was found
if AMIRAMESH
ALLTESTS += amirameshtest
endif

if HAVE_ALBERTA
ALLTESTS += gmshtest-alberta
endif

# Currently, Star-CD files can only be read into UGGrid objects
if UG
ALLTESTS += starcdreadertest
endif

# whats that for?
AM_CPPFLAGS += -g -Wall -O0

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS)

vtktest_SOURCES = vtktest.cc
vtktest_LDADD = $(LOCAL_LIBS)

vtksequencetest_SOURCES = vtksequencetest.cc
vtksequencetest_LDADD = $(LOCAL_LIBS)

gnuplottest_SOURCES = gnuplottest.cc
gnuplottest_LDADD = $(LOCAL_LIBS)

subsamplingvtktest_SOURCES = subsamplingvtktest.cc test-linking.cc
subsamplingvtktest_LDADD = $(LOCAL_LIBS)

amirameshtest_SOURCES = amirameshtest.cc
amirameshtest_CXXFLAGS = $(AMIRAMESH_CPPFLAGS) $(UG_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(PSURFACE_CPPFLAGS)
amirameshtest_LDADD = $(LOCAL_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)

gmshtest_SOURCES = gmshtest.cc
gmshtest_CXXFLAGS = $(UG_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
gmshtest_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

gmshtest_alberta_SOURCES = gmshtest-alberta.cc
gmshtest_alberta_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS)
gmshtest_alberta_LDADD = $(LOCAL_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(GRAPE_LDFLAGS) $(GRAPE_LIBS)

gmshtest_all_SOURCES = gmshtest-alberta.cc
gmshtest_all_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
gmshtest_all_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

starcdreadertest_SOURCES = starcdreadertest.cc
starcdreadertest_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)
starcdreadertest_LDADD = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp *.data sgrid*.am *.pvtu *.pvtp *.pvd
