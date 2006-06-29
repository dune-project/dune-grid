# $Id$

if BUILD_DOCS
 PAGES =
endif

SUBDIRS = refelements

all: $(PAGES)

# install the pages in the distribution
docdir=$(datadir)/doc/dune-grid/appl
doc_DATA = $(PAGES)

# place in website
BASEDIR=../..
CURDIR=doc/appl

include $(top_srcdir)/am/webstuff

CLEANFILES = $(PAGES)

include $(top_srcdir)/am/global-rules
