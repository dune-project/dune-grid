ALLTESTS = vtktest gnuplottest vtksequencetest subsamplingvtktest gmshtest

GRIDDIM=2
GRIDTYPE=YASPGRID

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

# run amiramesh test only if amiramesh was found
if AMIRAMESH
ALLTESTS += amirameshtest
endif

if ALBERTA
ALLTESTS += gmshtest-alberta2d
ALLTESTS += gmshtest-alberta3d
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

ALLTESTS += conformvolumevtktest
conformvolumevtktest_SOURCES = conformvolumevtktest.cc

ALLTESTS += nonconformboundaryvtktest
nonconformboundaryvtktest_SOURCES = nonconformboundaryvtktest.cc

vtktest_SOURCES = vtktest.cc

vtksequencetest_SOURCES = vtksequencetest.cc

gnuplottest_SOURCES = gnuplottest.cc

subsamplingvtktest_SOURCES = subsamplingvtktest.cc test-linking.cc
subsamplingvtktest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)
subsamplingvtktest_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)
subsamplingvtktest_LDADD =			\
	$(DUNEMPILIBS)				\
	$(LDADD)

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

gmshtest_SOURCES  = gmshtest.cc
gmshtest_CPPFLAGS = $(AM_CPPFLAGS)  $(ALL_PKG_CPPFLAGS)
gmshtest_LDFLAGS  = $(AM_LDFLAGS)   $(ALL_PKG_LDFLAGS)
gmshtest_LDADD    = $(ALL_PKG_LIBS) $(LDADD)

gmshtest_alberta2d_SOURCES  = gmshtest.cc
gmshtest_alberta2d_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS) -DGRIDDIM=$(ALBERTA_DIM)
gmshtest_alberta2d_LDFLAGS  = $(AM_LDFLAGS)  $(ALBERTA_LDFLAGS)  $(GRAPE_LDFLAGS)
gmshtest_alberta2d_LDADD    =	$(GRAPE_LIBS)	 $(ALBERTA_LIBS)		 $(LDADD)

gmshtest_alberta3d_SOURCES  = gmshtest.cc
gmshtest_alberta3d_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS) $(GRAPE_CPPFLAGS)
gmshtest_alberta3d_LDFLAGS  = $(AM_LDFLAGS)  $(ALBERTA3D_LDFLAGS)  $(GRAPE_LDFLAGS)
gmshtest_alberta3d_LDADD    =	$(GRAPE_LIBS)	 $(ALBERTA3D_LIBS)		 $(LDADD)

gmshtest_alugrid_SOURCES  = gmshtest.cc
gmshtest_alugrid_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(GRAPE_CPPFLAGS) $(GRIDDIM_CPPFLAGS)
gmshtest_alugrid_LDFLAGS  = $(AM_LDFLAGS)	 $(ALUGRID_LDFLAGS)	 $(GRAPE_LDFLAGS)
gmshtest_alugrid_LDADD    =	$(GRAPE_LIBS)	 $(ALUGRID_LIBS)		 $(LDADD)

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
