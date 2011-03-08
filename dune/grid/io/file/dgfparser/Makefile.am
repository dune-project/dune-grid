SUBDIRS = blocks test utils

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = dgfparser.cc

libdgfparser_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(DUNE_CPPFLAGS)
libdgfparser_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(DUNE_LDFLAGS)
libdgfparser_la_LIBADD = blocks/libdgfparserblocks.la \
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)

dgfparserdir = $(includedir)/dune/grid/io/file/dgfparser
dgfparser_HEADERS = dgfalberta.hh  dgfalu.cc  dgfexception.hh  \
		    dgfug.cc   dgfalu.hh  dgfug.hh \
		    dgfparser.hh  dgfgeogrid.hh \
		    dgfwriter.hh  dgfyasp.hh  dgfgridtype.hh \
		    entitykey.hh  dgfs.hh  entitykey_inline.hh  \
		    dgfoned.cc  dgfoned.hh

include $(top_srcdir)/am/global-rules
