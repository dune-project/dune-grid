AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(MPI_LDFLAGS) $(MPI_LIBS)
AM_CPPFLAGS =  @AM_CPPFLAGS@ $(MPI_CPPFLAGS) -DSRCDIR="\"$(srcdir)\""

if ALUGRID
  TESTALU = testalu
endif

if ALBERTA
  TESTALBERTA = testalberta
endif

if UG
  TESTUG = testug test2ug
endif

if GRAPE
  VIEWPROGS = viewdgf
endif

ALLTESTS = $(TESTALU) $(TESTALBERTA) testsgrid testyasp testoned $(TESTUG)

# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS) tester $(VIEWPROGS)

# list of tests to run
TESTS = $(ALLTESTS)

tester_SOURCES = main.cc
tester_CXXFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
tester_LDADD = $(LOCAL_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(MPI_LDFLAGS) $(UG_LIBS) $(UG_LDFLAGS)

viewdgf_SOURCES = viewdgf.cc
viewdgf_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
viewdgf_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

if ALUGRID
testalu_CXXFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(GRAPE_CPPFLAGS)
testalu_SOURCES = testalu.cc
testalu_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)
endif

if ALBERTA
testalberta_CXXFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS)
testalberta_SOURCES = testalberta.cc
testalberta_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS)
endif

testsgrid_CXXFLAGS = $(AM_CPPFLAGS)
testsgrid_SOURCES = testsg.cc
testsgrid_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS)

testyasp_CXXFLAGS = $(AM_CPPFLAGS)
testyasp_SOURCES = testyasp.cc
testyasp_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS)

testoned_CXXFLAGS = $(AM_CPPFLAGS)
testoned_SOURCES = testoned.cc
testoned_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS)

if UG
testug_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)
testug_SOURCES = testug.cc
testug_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS) $(UG_LDFLAGS) $(UG_LIBS)
test2ug_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)
test2ug_SOURCES = test2ug.cc
test2ug_LDADD = $(LOCAL_LIBS) $(GRAPE_LIBS) $(GRAPE_LDFLAGS) $(UG_LDFLAGS) $(UG_LIBS)
endif

GRIDFILES = alberta-testgrid-3-3.al \
  simplex-testgrid-2-2.dgf simplex-testgrid-3-3-large.dgf simplex-testgrid-3-3.dgf \
  cube-testgrid-3-3-large.dgf cube-testgrid-3-3.dgf

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = examplegrid10a.dgf examplegrid10b.dgf examplegrid10.dgf		\
	examplegrid11.dgf examplegrid12.dgf examplegrid1c.dgf					\
	examplegrid1gen.dgf examplegrid1s.dgf examplegrid2a.dgf					\
	examplegrid2b.dgf examplegrid2c.dgf examplegrid2d.dgf					\
	examplegrid2e.dgf examplegrid5.dgf examplegrid6.dgf examplegrid7.dgf	\
	examplegrid9.dgf

CLEANFILES = examplegrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules

