# $Id$

if UG

noinst_LTLIBRARIES = libuggrid.la

libuggrid_la_SOURCES  = uggrid.cc uggridentity.cc boundaryextractor.cc boundaryextractor.hh \
                        uggridindexsets.cc ugintersectionit.cc uggridfactory.cc

libuggrid_la_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(MPI_CPPFLAGS)
libuggrid_la_LIBADD   = $(UG_LDFLAGS) $(UG_LIBS) $(DUNE_COMMON_LIBS) $(MPI_LDFLAGS) $(MPI_LIBS)

endif

headercheck_IGNORE = uggridentity.hh uggridentitypointer.hh uggridgeometry.hh \
                     uggridhieriterator.hh uggridindexsets.hh uggridleafiterator.hh \
                     uggridleveliterator.hh ugincludes.hh ugintersectionit.hh \
                     ugwrapper.hh

uggriddir = $(includedir)/dune/grid/uggrid/
uggrid_HEADERS = uggridfactory.hh uggridgeometry.cc \
  uggridentitypointer.hh uggridhieriterator.cc \
  ugwrapper.hh uggridentity.hh uggridgeometry.hh \
  uggridhieriterator.hh uggridleveliterator.hh ugincludes.hh \
  ugintersectionit.hh uggridindexsets.hh \
  uggridleafiterator.hh uggridrenumberer.hh \
  ug_undefs.hh ug_undefs_lgm_seq.hh

# tricks like undefAllMacros.pl don't have to be shipped, have they?
# Oliver Sander: It depends. undefAllMacros.pl creates the file ug_undefs.hh
# which contains an #undef line for each macro encountered in the UG
# headers which DUNE includes.  This should be redone every time UG
# changes, and currently I do it and then check in the new ug_undefs.hh.
# An alternative would be to have the DUNE-build system call undefAllMacros
# for the UG it is configured with.  Then it needed to be shipped.


include $(top_srcdir)/am/global-rules
