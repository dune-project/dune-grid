# $Id$

if PSURFACE

noinst_LTLIBRARIES = libamiramesh.la

libamiramesh_la_SOURCES = psurface.cc

libamiramesh_la_CPPFLAGS = $(AM_CPPFLAGS) \
                           $(PSURFACE_CPPFLAGS) \
                           $(AMIRAMESH_CPPFLAGS)
libamiramesh_la_LDFLAGS  = $(AM_LDFLAGS)
libamiramesh_la_LIBADD   = $(PSURFACE_LIBS) \
                           $(AMIRAMESH_LIBS)

endif

amirameshdir = $(includedir)/dune/grid/io/file/amiramesh
amiramesh_HEADERS = amirameshreader.cc amirameshwriter.cc \
                    psurface.h

include $(top_srcdir)/am/global-rules
