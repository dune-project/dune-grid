SUBDIRS = blocks test utils

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = dgfparser.cc dgfug.cc

libdgfparser_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(DUNE_CPPFLAGS)			\
	$(UG_CPPFLAGS)
libdgfparser_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(DUNE_LDFLAGS)				\
	$(UG_LDFLAGS)
libdgfparser_la_LIBADD = blocks/libdgfparserblocks.la \
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)				\
	$(UG_LIBS)

dgfparserdir = $(includedir)/dune/grid/io/file/dgfparser
dgfparser_HEADERS = dgfalberta.hh  dgfalu.cc  dgfexception.hh  \
		    dgfalu.hh  dgfug.hh \
		    dgfparser.hh  dgfgeogrid.hh \
		    dgfwriter.hh  dgfyasp.hh \
		    entitykey.hh  dgfs.hh  entitykey_inline.hh  \
		    dgfoned.hh dgfgridfactory.hh \
		    macrogrid.hh  gridptr.hh  parser.hh

include $(top_srcdir)/am/global-rules
