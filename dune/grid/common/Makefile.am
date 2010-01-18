# $Id$

noinst_LTLIBRARIES = libgridcommon.la
libgridcommon_la_SOURCES = referenceelements.cc
libgridcommon_la_LIBADD = quadraturerules/libquadraturerules.la \
  $(DUNE_COMMON_LIBS)

commondir = $(includedir)/dune/grid/common
common_HEADERS = virtualrefinement.cc capabilities.hh \
  entity.hh entitypointer.hh geometry.hh grid.hh \
  hierarchiciterator.hh indexstack.hh intersection.hh intersectioniterator.hh \
  leafiterator.hh leveliterator.hh refinement.hh \
  virtualrefinement.hh referenceelements.hh genericreferenceelements.hh \
  boundarysegment.hh gridinfo.hh indexidset.hh mapper.hh \
  mcmgmapper.hh scsgmapper.hh universalmapper.hh sizecache.hh \
  intersectioniteratorwrapper.hh quadraturerules.hh \
  datahandleif.hh gridfactory.hh \
  gridview.hh defaultgridview.hh \
  gridenums.hh adaptcallback.hh \
  dynamicsubindexid.hh boundaryprojection.hh \
  topologyfactory.hh

SUBDIRS = refinement test quadraturerules

dist_noinst_DATA = cube-to-tet-6.png

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/inkscape.am
