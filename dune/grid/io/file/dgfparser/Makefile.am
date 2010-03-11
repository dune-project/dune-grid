SUBDIRS = test utils

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = dgfparserblocks.cc dgfprojectionblock.cc dgfparser.cc
libdgfparser_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(DUNE_CPPFLAGS)
libdgfparser_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(DUNE_LDFLAGS)
libdgfparser_la_LIBADD =			\
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)

dgfparserdir = $(includedir)/dune/grid/io/file/dgfparser
dgfparser_HEADERS = dgfwriter.hh \
dgfalberta.hh  dgfoned.hh    dgfparserblocks.hh  dgfug.hh    entitykey.hh entitykey_inline.hh \
dgfalu.hh      dgfparser.hh  dgfyasp.hh  dgfs.hh \
dgfalu.cc      dgfoned.cc    dgfug.cc \
dgfgridtype.hh  dgfexception.hh  dgfgeogrid.hh  dgfprojectionblock.hh

include $(top_srcdir)/am/global-rules
