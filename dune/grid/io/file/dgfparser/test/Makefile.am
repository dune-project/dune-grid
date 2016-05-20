AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS) -DHAVE_DUNE_GRID=1
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)
LDADD = 					\
	$(DUNEMPILIBS)				\
	@LDADD@

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

if ALBERTA
  TESTALBERTA = testalberta
endif

if UG
  TESTUG = testug
endif

ALLTESTS = $(TESTALBERTA) testyasp testdgfyaspoffset testoned $(TESTUG)

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)
EXTRA_PROGRAMS = tester

# list of tests to run
TESTS = $(ALLTESTS)

tester_SOURCES = test-dgf.cc
tester_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
tester_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALL_PKG_LDFLAGS)
tester_LDADD =					\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

if ALBERTA
testalberta_SOURCES = test-dgf.cc
testalberta_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA2D_CPPFLAGS)			\
	-DALBERTAGRID -DGRIDDIM=2 -DWORLDDIM=2
testalberta_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA2D_LDFLAGS)
testalberta_LDADD =				\
	$(ALBERTA2D_LIBS)	\
	$(LDADD)
endif

testyasp_SOURCES = test-dgf.cc
testyasp_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DYASPGRID -DGRIDDIM=3

testdgfyaspoffset_SOURCES = test-dgf.cc
testdgfyaspoffset_CPPFLAGS = $(AM_CPPFLAGS)	\
	-DYASPGRID -DGRIDDIM=3			\
	-DCOMPLETE_GRID_TYPE="YaspGrid<2, EquidistantOffsetCoordinates<double, 2> >"

testoned_SOURCES = test-dgf.cc
testoned_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DONEDGRID -DGRIDDIM=1

if UG
testug_SOURCES = test-dgf.cc
testug_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(UG_CPPFLAGS)				\
	-DUGGRID -DGRIDDIM=3
testug_LDFLAGS = $(AM_LDFLAGS)			\
	$(UG_LDFLAGS)
testug_LDADD =					\
	$(UG_LIBS)				\
	$(LDADD)
endif

CLEANFILES = dgfparser.log

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
