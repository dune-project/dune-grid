# $Id$

noinst_LTLIBRARIES = libgridcommon.la
# construct only from sub-libraries
libgridcommon_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libgridcommon_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libgridcommon_la_LIBADD = quadraturerules/libquadraturerules.la \
  $(DUNE_LIBS)

commondir = $(includedir)/dune/grid/common
common_HEADERS = virtualrefinement.cc capabilities.hh \
  entity.hh entitypointer.hh geometry.hh grid.hh \
  hierarchiciterator.hh indexstack.hh intersection.hh intersectioniterator.hh \
  leafiterator.hh leveliterator.hh entityiterator.hh refinement.hh \
  virtualrefinement.hh genericreferenceelements.hh \
  boundarysegment.hh gridinfo.hh indexidset.hh mapper.hh \
  mcmgmapper.hh scsgmapper.hh universalmapper.hh sizecache.hh \
  intersectioniteratorwrapper.hh quadraturerules.hh \
  datahandleif.hh gridfactory.hh \
  gridview.hh defaultgridview.hh \
  gridenums.hh adaptcallback.hh \
  dynamicsubindexid.hh boundaryprojection.hh \
  topologyfactory.hh exceptions.hh

SUBDIRS = refinement test quadraturerules

dist_noinst_DATA = cube-to-tet-6.png

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/inkscape.am
