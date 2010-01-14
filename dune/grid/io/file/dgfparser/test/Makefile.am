AM_LDFLAGS = @AM_LDFLAGS@			\
	$(ALL_PKG_LDFLAGS)			\
	$(DUNEMPILDFLAGS)
LDADD = $(DUNEMPILIBS)
AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS)			\
	-DSRCDIR="\"$(srcdir)\""

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
check_PROGRAMS = $(ALLTESTS) $(VIEWPROGS)
EXTRA_PROGRAMS = tester

# list of tests to run
TESTS = $(ALLTESTS)

tester_SOURCES = main.cc
tester_CXXFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
tester_LDADD = $(LOCAL_LIBS)			\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)	\
	$(GRAPE_LDFLAGS) $(GRAPE_LIBS)		\
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)	\
	$(UG_LDFLAGS) $(UG_LIBS)		\
	$(LDADD)

viewdgf_SOURCES = viewdgf.cc
viewdgf_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
viewdgf_LDADD = $(LOCAL_LIBS)			\
	$(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)	\
	$(LDADD)

if ALUGRID
testalu_CXXFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) -DALUGRID_SIMPLEX -DGRIDDIM=3
testalu_SOURCES = main.cc
testalu_LDADD = $(LOCAL_LIBS)			\
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)	\
	$(LDADD)
endif

if ALBERTA
testalberta_CXXFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) -DALBERTAGRID -DGRIDDIM=2
testalberta_SOURCES = main.cc
testalberta_LDADD = $(LOCAL_LIBS)		\
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)	\
	$(LDADD)
endif

testsgrid_CXXFLAGS = $(AM_CPPFLAGS) -DSGRID -DGRIDDIM=3
testsgrid_SOURCES = main.cc
testsgrid_LDADD = $(LOCAL_LIBS) $(LDADD)

testyasp_CXXFLAGS = $(AM_CPPFLAGS) -DYASPGRID -DGRIDDIM=3
testyasp_SOURCES = main.cc
testyasp_LDADD = $(LOCAL_LIBS) $(LDADD)

testoned_CXXFLAGS = $(AM_CPPFLAGS) -DONEDGRID -DGRIDDIM=1
testoned_SOURCES = main.cc
testoned_LDADD = $(LOCAL_LIBS) $(LDADD)

if UG
testug_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) -DUGGRID -DGRIDDIM=3
testug_SOURCES = main.cc
testug_LDADD = $(LOCAL_LIBS)			\
	$(UG_LDFLAGS) $(UG_LIBS)		\
	$(LDADD)
endif

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = examplegrid10a.dgf examplegrid10b.dgf examplegrid10.dgf		\
	examplegrid11.dgf examplegrid12.dgf examplegrid1c.dgf					\
	examplegrid1gen.dgf examplegrid1s.dgf examplegrid2a.dgf					\
	examplegrid2b.dgf examplegrid2c.dgf examplegrid2d.dgf					\
	examplegrid2e.dgf examplegrid5.dgf examplegrid6.dgf examplegrid7.dgf	\
	examplegrid9.dgf test3d.dgf test2d.dgf test1d.dgf

CLEANFILES = examplegrid*.dgf.* dgfparser.log test?d.dgf.*

include $(top_srcdir)/am/global-rules
