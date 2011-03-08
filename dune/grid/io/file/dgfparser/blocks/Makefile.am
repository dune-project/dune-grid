noinst_LTLIBRARIES = libdgfparserblocks.la

libdgfparserblocks_la_SOURCES = basic.cc  boundarydom.cc  boundaryseg.cc \
				cube.cc  dim.cc  general.cc  gridparameter.cc \
				interval.cc  periodicfacetrans.cc  projection.cc \
				simplex.cc  simplexgeneration.cc  vertex.cc

libdgfparserblocks_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)			\
	$(DUNE_CPPFLAGS)
libdgfparserblocks_la_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(DUNE_LDFLAGS)
libdgfparserblocks_la_LIBADD =			\
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)

blocksdir = $(includedir)/dune/grid/io/file/dgfparser/blocks
blocks_HEADERS = basic.hh  boundarydom.hh  boundaryseg.hh \
		 cube.hh  dim.hh  general.hh  gridparameter.hh \
		 interval.hh  periodicfacetrans.hh  projection.hh \
		 simplex.hh  simplexgeneration.hh  vertex.hh

include $(top_srcdir)/am/global-rules
