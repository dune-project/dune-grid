# $Id$

SUBDIRS = alberta amiramesh test dgfparser

iofiledir = $(includedir)/dune/grid/io/file
iofile_HEADERS = amirameshreader.hh amirameshwriter.hh \
 albertawriter.hh alugridwriter.hh

noinst_LTLIBRARIES = libfileio.la

if AMIRAMESH
if UG
  AMIRALIB = amiramesh/libugamiramesh.la
endif
endif

libfileio_la_SOURCES =
libfileio_la_LIBADD = $(AMIRALIB) -lcommon

include $(top_srcdir)/am/global-rules
