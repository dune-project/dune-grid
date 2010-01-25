AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS)			\
	-DSRCDIR="\"$(srcdir)\""
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)
LDADD = 					\
	$(DUNEMPILIBS)				\
	@LDADD@

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
	$(ALL_PKG_CPPFLAGS)
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
	-DALUGRID_SIMPLEX -DGRIDDIM=3
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

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = examplegrid10a.dgf examplegrid10b.dgf examplegrid10.dgf		\
	examplegrid11.dgf examplegrid12.dgf examplegrid1c.dgf					\
	examplegrid1gen.dgf examplegrid1s.dgf examplegrid2a.dgf					\
	examplegrid2b.dgf examplegrid2c.dgf examplegrid2d.dgf					\
	examplegrid2e.dgf examplegrid5.dgf examplegrid6.dgf examplegrid7.dgf	\
	examplegrid9.dgf test3d.dgf test2d.dgf test1d.dgf \
	example-projection.dgf

CLEANFILES = examplegrid*.dgf.* dgfparser.log test?d.dgf.*

include $(top_srcdir)/am/global-rules
