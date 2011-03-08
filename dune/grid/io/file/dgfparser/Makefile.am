SUBDIRS = test utils

noinst_LTLIBRARIES = libdgfparser.la

libdgfparser_la_SOURCES = dgfbasicblock.cc  dgfcubeblock.cc  \
			  dgfgridparameterblock.cc  dgfprojectionblock.cc  \
			  dgfvertexblock.cc  dgfboundarydomblock.cc \
			  dgfdimblock.cc  dgfintervalblock.cc  dgfsimplexblock.cc \
			  dgfboundarysegblock.cc  dgfgeneralblock.cc  \
			  dgfperiodicfacetransblock.cc  dgfsimplexgenerationblock.cc \
			  dgfparser.cc

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
dgfparser_HEADERS = dgfalberta.hh  dgfdimblock.hh  dgfsimplexgenerationblock.hh \
		    dgfalu.cc  dgfexception.hh  dgfug.cc   dgfalu.hh \
		    dgfug.hh  dgfgeneralblock.hh \
		    dgfparser.hh  dgfbasicblock.hh  dgfgeogrid.hh  \
		    dgfvertexblock.hh \
		    dgfperiodicfacetransblock.hh  dgfwriter.hh \
		    dgfboundarydomblock.hh  dgfgridparameterblock.hh \
		    dgfyasp.hh  dgfgridtype.hh  dgfprojectionblock.hh \
		    entitykey.hh  dgfboundarysegblock.hh  dgfs.hh \
		    entitykey_inline.hh  dgfintervalblock.hh \
		    dgfcubeblock.hh  dgfoned.cc  dgfsimplexblock.hh  dgfoned.hh

include $(top_srcdir)/am/global-rules
