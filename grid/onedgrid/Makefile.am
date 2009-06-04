# $Id$

noinst_LTLIBRARIES = libonedgrid.la

libonedgrid_la_SOURCES = onedgrid.cc onedgridfactory.cc nulliteratorfactory.cc
libonedgrid_la_LIBADD = $(DUNE_COMMON_LIBS)

onedgriddir = $(includedir)/dune/grid/onedgrid/
onedgrid_HEADERS = nulliteratorfactory.hh  onedgridentity.hh \
   onedgridentitypointer.hh onedgridfactory.hh onedgridgeometry.hh  onedgridhieriterator.hh \
   onedgridindexsets.hh  onedgridleafiterator.hh  onedgridleveliterator.hh \
   onedgridlist.hh  onedintersectionit.hh

include $(top_srcdir)/am/global-rules
