# $Id$

yaspgriddir = $(includedir)/dune/grid/yaspgrid/
yaspgrid_HEADERS = coordinates.hh \
                   factory.hh \
                   sequences.hh \
                   torus.hh \
                   yaspgridentity.hh \
                   yaspgridentityseed.hh \
                   yaspgridentitypointer.hh \
                   yaspgridgeometry.hh \
                   yaspgridhierarchiciterator.hh \
                   yaspgrididset.hh \
                   yaspgridindexsets.hh \
                   yaspgridintersection.hh \
                   yaspgridintersectioniterator.hh \
                   yaspgridleveliterator.hh \
                   ygrid.hh

# The header yaspgrid.hh declares a few global variables.  These are used
# in most other headers, and therefore those cannot currently pass the headercheck.
headercheck_IGNORE = yaspgridentity.hh \
                     yaspgridentityseed.hh \
                     yaspgridentitypointer.hh \
                     yaspgridgeometry.hh \
                     yaspgridhierarchiciterator.hh \
                     yaspgrididset.hh \
                     yaspgridindexsets.hh \
                     yaspgridintersection.hh \
                     yaspgridintersectioniterator.hh \
                     yaspgridleveliterator.hh

EXTRA_DIST = CMakeLists.txt grid.fig grid.eps grid.png subgrid.fig subgrid.eps subgrid.png

include $(top_srcdir)/am/global-rules
