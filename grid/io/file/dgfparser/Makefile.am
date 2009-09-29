SUBDIRS = test utils

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = dgfparserblocks.cc dgfparser.cc
libdgfparser_la_LIBADD = $(DUNE_COMMON_LIBS)

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS)

AM_LDFLAGS=

dgfparserdir = $(includedir)/dune/grid/io/file/dgfparser
dgfparser_HEADERS = dgfwriter.hh \
dgfalberta.hh  dgfoned.hh    dgfparserblocks.hh  dgfug.hh    entitykey.hh entitykey_inline.hh \
dgfalu.hh      dgfparser.hh  dgfyasp.hh  dgfs.hh gridtype.hh \
dgfalu.cc      dgfoned.cc    dgfug.cc   dgfyasp.cc   dgfs.cc \
dgfgridtype.hh  dgfexception.hh  dgfgeogrid.hh

include $(top_srcdir)/am/global-rules
