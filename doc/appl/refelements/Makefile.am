# $Id$

# include this into website
EXTRAINSTALL=hexahedron.png line.png prism.png pyramid.png \
  quadrilateral.png tetrahedron.png triangle.png\
  hexahedron_edges.png tetrahedron_edges.png prism_edges.png pyramid_edges.png

# place in website
BASEDIR=../../..
CURDIR=doc/appl/refelements

# install docs
docdir=$(datadir)/doc/dune-grid/appl/refelements
if BUILD_DOCS
DOCFILES = $(EXTRAINSTALL)
EXTRA_DIST= $(EXTRAINSTALL)
endif

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
