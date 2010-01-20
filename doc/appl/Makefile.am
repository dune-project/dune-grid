# $Id$

if BUILD_DOCS
PAGES =
endif

SUBDIRS = refelements

# place in website
BASEDIR=../..
CURDIR=doc/appl

# install the html pages
docdir=$(datadir)/doc/dune-grid/appl
if BUILD_DOCS
DOCFILES = $(PAGES)
EXTRA_DIST = $(PAGES)
endif

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
