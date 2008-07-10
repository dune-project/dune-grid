SUBDIRS = test

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = entitykey.cc dgfparserblocks.cc dgfparser.cc
libdgfparser_la_LIBADD = $(DUNE_COMMON_LIBS)

AM_LDFLAGS=

dgfparserdir = $(includedir)/dune/grid/io/file/dgfparser
dgfparser_HEADERS = \
dgfalberta.hh  dgfoned.hh    dgfparserblocks.hh  dgfug.hh    entitykey.hh \
dgfalu.hh      dgfparser.hh  dgfyasp.hh  dgfs.hh gridtype.hh \
dgfalberta.cc  dgfalu.cc  dgfoned.cc  dgfparser.cc  dgfug.cc dgfyasp.cc dgfs.cc \
dgfgridtype.hh  dgfexception.hh

include $(top_srcdir)/am/global-rules
