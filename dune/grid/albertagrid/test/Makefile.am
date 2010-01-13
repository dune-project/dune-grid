EXTRA_DIST =
check_PROGRAMS =

TESTS =
XFAIL_TESTS =

GRIDDIM=2

EXTRA_DIST += svgout.cc

if ALBERTA
check_PROGRAMS += svgout
svgout_SOURCES = svgout.cc
svgout_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
svgout_LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# test program will return the special exit code 77 ("skipped") if alberta is
# not available
TESTS += test-alberta3d-refine
XFAIL_TESTS += test-alberta3d-refine
check_PROGRAMS += test-alberta3d-refine
test_alberta3d_refine_SOURCES = test-alberta3d-refine.cc
test_alberta3d_refine_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
test_alberta3d_refine_LDADD = $(LOCAL_LIBS) \
	$(ALBERTA_LDFLAGS) -lalberta_3d $(ALBERTA_BASE_LIBS)
endif

include $(top_srcdir)/am/global-rules
