# $Id$

if INKSCAPE
INKSCAPE_SUBDIRS = refinement
endif INKSCAPE
SUBDIRS = appl $(INKSCAPE_SUBDIRS) doxygen

# setting like in dune-web
BASEDIR=..
CURDIR=doc

if ! BUILD_DOCS
# add tag to notify that dist has been build without documentation
dist-hook:
	echo "# No documentation included in distribution! " > $(distdir)/$(DOCUMENTATION_TAG_FILE)
endif

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
