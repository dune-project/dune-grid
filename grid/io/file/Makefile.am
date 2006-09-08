# $Id$

SUBDIRS = amiramesh test dgfparser vtk

iofiledir = $(includedir)/dune/grid/io/file
iofile_HEADERS = amirameshreader.hh amirameshwriter.hh

noinst_LTLIBRARIES = libfileio.la

if AMIRAMESH
if UG
  AMIRALIB = amiramesh/libugamiramesh.la
endif
endif

libfileio_la_SOURCES =
libfileio_la_LIBADD = $(AMIRALIB) $(DUNE_COMMON_LIBS)

include $(top_srcdir)/am/global-rules
