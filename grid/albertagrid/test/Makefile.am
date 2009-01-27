if ALBERTA
noinst_PROGRAMS = svgout
endif

GRIDDIM=2

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
LDADD = $(LOCAL_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

svgout_SOURCES = svgout.cc

EXTRA_DIST = svgout.cc

include $(top_srcdir)/am/global-rules
