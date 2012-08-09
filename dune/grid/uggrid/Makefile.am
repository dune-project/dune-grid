# $Id$

if UG

noinst_LTLIBRARIES = libuggrid.la

libuggrid_la_SOURCES  = uggridintersections.cc uggrid.cc uggridentity.cc boundaryextractor.cc boundaryextractor.hh \
                        uggridindexsets.cc uggridfactory.cc uggridgeometry.cc uggridhieriterator.cc

libuggrid_la_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(UG_CPPFLAGS)				\
	$(DUNEMPICPPFLAGS)
libuggrid_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(UG_LDFLAGS)				\
	$(DUNEMPILDFLAGS)
libuggrid_la_LIBADD =				\
	$(UG_LIBS)				\
	$(DUNEMPILIBS)				\
	$(DUNE_LIBS)

endif

headercheck_IGNORE = uggridentity.hh uggridentitypointer.hh uggridentityseed.hh uggridfactory.hh uggridgeometry.hh \
                     uggridhieriterator.hh uggridindexsets.hh uggridleafiterator.hh \
                     uggridleveliterator.hh uggridlocalgeometry.hh uggridrenumberer.hh \
                     ugincludes.hh uggridintersections.hh \
                     uggridintersectioniterators.hh ugwrapper.hh

uggriddir = $(includedir)/dune/grid/uggrid/
uggrid_HEADERS = uggridfactory.hh uggridentitypointer.hh \
  uggridentityseed.hh uggridentity.hh uggridgeometry.hh \
  uggridlocalgeometry.hh \
  uggridhieriterator.hh uggridleveliterator.hh ugincludes.hh \
  uggridintersections.hh uggridintersectioniterators.hh uggridindexsets.hh \
  uggridleafiterator.hh uggridrenumberer.hh \
  ug_undefs.hh ugwrapper.hh

# tricks like undefAllMacros.pl don't have to be shipped, have they?
# Oliver Sander: It depends. undefAllMacros.pl creates the file ug_undefs.hh
# which contains an #undef line for each macro encountered in the UG
# headers which DUNE includes.  This should be redone every time UG
# changes, and currently I do it and then check in the new ug_undefs.hh.
# An alternative would be to have the DUNE-build system call undefAllMacros
# for the UG it is configured with.  Then it needed to be shipped.


include $(top_srcdir)/am/global-rules
