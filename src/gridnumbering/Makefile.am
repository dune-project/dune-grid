noinst_HEADERS =				\
	to-vtk.hh

bin_PROGRAMS =

AM_CPPFLAGS = @AM_CPPFLAGS@			\
	$(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@			\
	$(DUNEMPILDFLAGS)
LDADD =						\
	$(DUNEMPILIBS)				\
	@LDADD@


bin_PROGRAMS += gridnumbering-yasp-1d
gridnumbering_yasp_1d_SOURCES = yasp.cc
gridnumbering_yasp_1d_CPPFLAGS = $(AM_CPPFLAGS)	\
	-DDIMENSION=1

bin_PROGRAMS += gridnumbering-yasp-2d
gridnumbering_yasp_2d_SOURCES = yasp.cc
gridnumbering_yasp_2d_CPPFLAGS = $(AM_CPPFLAGS)	\
	-DDIMENSION=2

bin_PROGRAMS += gridnumbering-yasp-3d
gridnumbering_yasp_3d_SOURCES = yasp.cc
gridnumbering_yasp_3d_CPPFLAGS = $(AM_CPPFLAGS)	\
	-DDIMENSION=3

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
