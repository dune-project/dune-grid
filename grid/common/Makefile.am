# $Id$

noinst_LTLIBRARIES = libgridcommon.la
libgridcommon_la_SOURCES = referenceelements.cc quadraturerules.cc
libgridcommon_la_LIBADD = -lcommon

commondir = $(includedir)/dune/grid/common
common_HEADERS = virtualrefinement.cc capabilities.hh \
  defaultindexsets.hh entity.hh entitypointer.hh geometry.hh grid.hh \
  hierarchiciterator.hh indexstack.hh intersectioniterator.hh \
  leafiterator.hh leveliterator.hh refinement.hh \
  virtualrefinement.hh referenceelements.hh gridpart.hh \
  boundarysegment.hh gridinfo.hh indexidset.hh mapper.hh \
  mcmgmapper.hh scsgmapper.hh universalmapper.hh sizecache.hh \
  intersectioniteratorwrapper.hh quadraturerules.hh \
  datahandleif.hh

SUBDIRS = refinement test

include $(top_srcdir)/am/global-rules
