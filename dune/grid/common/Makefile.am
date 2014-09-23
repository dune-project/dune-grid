# $Id$

commondir = $(includedir)/dune/grid/common
common_HEADERS = \
	adaptcallback.hh \
	backuprestore.hh \
	boundaryprojection.hh \
	boundarysegment.hh \
	capabilities.hh \
	datahandleif.hh \
	defaultgridview.hh \
	entity.hh \
	entityiterator.hh \
	entitypointer.hh \
	entityrangegenerators.hh \
	entityseed.hh \
	exceptions.hh \
	geometry.hh \
	grid.hh \
	gridenums.hh \
	gridfactory.hh \
	gridinfo.hh \
	gridview.hh \
	indexidset.hh \
	intersection.hh \
	intersectioniterator.hh \
	mcmgmapper.hh \
	mapper.hh \
	partitionset.hh \
	sizecache.hh \
	scsgmapper.hh \
	universalmapper.hh

SUBDIRS = test

if BUILD_DOCS
dist_noinst_DATA = cube-to-tet-6.png
endif

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/inkscape.am

EXTRA_DIST = CMakeLists.txt
