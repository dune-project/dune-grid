EXTRA_DIST =
check_PROGRAMS =

TESTS =
XFAIL_TESTS =

GRIDDIM=2

if ALBERTA
# test program will return the special exit code 77 ("skipped") if alberta is
# not available
TESTS += test-alberta3d-refine
# not a real test.
# it illustrates a problem with longest edge refinement
# that can be observed in certain alberta vervsion
EXTRA_PROGRAMS = test-alberta3d-refine
test_alberta3d_refine_SOURCES = test-alberta3d-refine.cc
test_alberta3d_refine_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(ALBERTA3D_CPPFLAGS)
test_alberta3d_refine_LDFLAGS = $(AM_LDFLAGS)	\
	$(ALBERTA3D_LDFLAGS)
test_alberta3d_refine_LDADD =			\
	$(ALBERTA3D_LIBS)			\
	$(LDADD)
endif

include $(top_srcdir)/am/global-rules
