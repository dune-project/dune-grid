ALLTESTS = vtktest gnuplottest vtksequencetest subsamplingvtktest gmshtest gmshtest-all

GRIDDIM=2
GRIDTYPE=YASPGRID

# run amiramesh test only if amiramesh was found
if AMIRAMESH
ALLTESTS += amirameshtest
endif

if ALBERTA
ALLTESTS += gmshtest-alberta
endif

if ALUGRID
ALLTESTS += gmshtest-alugrid
endif

# Currently, Star-CD files can only be read into UGGrid objects
if UG
ALLTESTS += starcdreadertest
endif

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS)

vtktest_SOURCES = vtktest.cc

vtksequencetest_SOURCES = vtksequencetest.cc

gnuplottest_SOURCES = gnuplottest.cc

subsamplingvtktest_SOURCES = subsamplingvtktest.cc test-linking.cc

amirameshtest_SOURCES = amirameshtest.cc
amirameshtest_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(AMIRAMESH_CPPFLAGS)			\
	$(UG_CPPFLAGS)				\
	$(ALBERTA_CPPFLAGS)			\
	$(ALUGRID_CPPFLAGS)			\
	$(PSURFACE_CPPFLAGS)
amirameshtest_LDFLAGS = $(AM_LDFLAGS)		\
	$(AMIRAMESH_LDFLAGS)			\
	$(UG_LDFLAGS)				\
	$(ALBERTA_LDFLAGS)			\
	$(ALUGRID_LDFLAGS)			\
	$(PSURFACE_LDFLAGS)
amirameshtest_LDADD =				\
	$(PSURFACE_LIBS)			\
	$(ALUGRID_LIBS)				\
	$(ALBERTA_LIBS)				\
	$(UG_LIBS)				\
	$(AMIRAMESH_LIBS)			\
	$(LDADD)

gmshtest_SOURCES = gmshtest.cc
gmshtest_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALL_PKG_CPPFLAGS)
gmshtest_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALL_PKG_LDFLAGS)
gmshtest_LDADD =				\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

gmshtest_alberta_SOURCES = gmshtest-alberta.cc
gmshtest_alberta_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS) -DGRIDDIM=$(ALBERTA_DIM)
gmshtest_alberta_CLDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA_LDFLAGS)			\
	$(GRAPE_LDFLAGS)
gmshtest_alberta_LDADD =			\
	$(GRAPE_LIBS)				\
	$(ALBERTA_LIBS)				\
	$(LDADD)

gmshtest_alugrid_SOURCES = gmshtest-alberta.cc
gmshtest_alugrid_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(GRAPE_CPPFLAGS) -DGRIDDIM=$(ALBERTA_DIM)
gmshtest_alugrid_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALUGRID_LDFLAGS)			\
	$(GRAPE_LDFLAGS)
gmshtest_alugrid_LDADD =			\
	$(GRAPE_LIBS)				\
	$(ALUGRID_LIBS)				\
	$(LDADD)

gmshtest_all_SOURCES = gmshtest-alberta.cc
gmshtest_all_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
gmshtest_all_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALL_PKG_LDFLAGS)
gmshtest_all_LDADD =				\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

starcdreadertest_SOURCES = starcdreadertest.cc
starcdreadertest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(UG_CPPFLAGS)
starcdreadertest_LDFLAGS = $(AM_LDFLAGS)	\
	$(UG_LDFLAGS)
starcdreadertest_LDADD =			\
	$(UG_LIBS)				\
	$(LDADD)

include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp *.data sgrid*.am *.pvtu *.pvtp *.pvd
