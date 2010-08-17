AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)
LDADD = 					\
	$(DUNEMPILIBS)				\
	@LDADD@

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

if ALUGRID
  TESTALU = testalu
endif

if ALBERTA
  TESTALBERTA = testalberta
endif

if UG
  TESTUG = testug
endif

if GRAPE
  VIEWPROGS = viewdgf
endif

ALLTESTS = $(TESTALU) $(TESTALBERTA) testsgrid testyasp testoned $(TESTUG)

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)
EXTRA_PROGRAMS = tester viewdgf

# list of tests to run
TESTS = $(ALLTESTS)

tester_SOURCES = main.cc
tester_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
tester_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALL_PKG_LDFLAGS)
tester_LDADD =					\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

viewdgf_SOURCES = viewdgf.cc
viewdgf_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
viewdgf_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALL_PKG_LDFLAGS)
viewdgf_LDADD =					\
	$(ALL_PKG_LIBS)				\
	$(LDADD)

if ALUGRID
testalu_SOURCES = main.cc
testalu_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALUGRID_CPPFLAGS)			\
	-DALUGRID_CUBE -DGRIDDIM=3
testalu_LDFLAGS = $(AM_LDFLAGS)			\
	$(ALUGRID_LDFLAGS)
testalu_LDADD =					\
	$(ALUGRID_LIBS)				\
	$(LDADD)
endif

if ALBERTA
testalberta_SOURCES = main.cc
testalberta_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA2D_CPPFLAGS)			\
	-DALBERTAGRID -DGRIDDIM=2 -DWORLDDIM=2
testalberta_LDFLAGS = $(AM_LDFLAGS)		\
	$(ALBERTA2D_LDFLAGS)
testalberta_LDADD =				\
	$(ALBERTA2D_LIBS)	\
	$(LDADD)
endif

testsgrid_SOURCES = main.cc
testsgrid_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DSGRID -DGRIDDIM=3

testyasp_SOURCES = main.cc
testyasp_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DYASPGRID -DGRIDDIM=3

testoned_SOURCES = main.cc
testoned_CPPFLAGS = $(AM_CPPFLAGS)		\
	-DONEDGRID -DGRIDDIM=1

if UG
testug_SOURCES = main.cc
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
