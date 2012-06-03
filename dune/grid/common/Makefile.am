# $Id$

commondir = $(includedir)/dune/grid/common
common_HEADERS = virtualrefinement.cc capabilities.hh \
  entity.hh entitypointer.hh geometry.hh grid.hh \
  hierarchiciterator.hh intersection.hh intersectioniterator.hh \
  leafiterator.hh leveliterator.hh entityiterator.hh refinement.hh \
  virtualrefinement.hh genericreferenceelements.hh \
  boundarysegment.hh gridinfo.hh indexidset.hh mapper.hh \
  mcmgmapper.hh scsgmapper.hh universalmapper.hh sizecache.hh \
  intersectioniteratorwrapper.hh quadraturerules.hh \
  datahandleif.hh gridfactory.hh \
  gridview.hh defaultgridview.hh \
  gridenums.hh adaptcallback.hh \
  dynamicsubindexid.hh boundaryprojection.hh \
  topologyfactory.hh exceptions.hh backuprestore.hh \
  geometryreference.hh

SUBDIRS = refinement test quadraturerules

dist_noinst_DATA = cube-to-tet-6.png

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/inkscape.am
