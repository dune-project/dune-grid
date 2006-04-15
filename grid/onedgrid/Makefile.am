# $Id$

noinst_LTLIBRARIES = libonedgrid.la

libonedgrid_la_SOURCES = onedgrid.cc
libonedgrid_la_CXXFLAGS = -DDUNE_EXPRESSIONTEMPLATES
libonedgrid_la_LIBADD = -lcommon

onedgriddir = $(includedir)/dune/grid/onedgrid/
onedgrid_HEADERS = \
   onedgridgeometry.hh onedgridentity.hh onedgridhieriterator.hh \
   onedgridleveliterator.hh onedintersectionit.hh onedgridentitypointer.hh \
   onedgridindexsets.hh onedgridleafiterator.hh

include $(top_srcdir)/am/global-rules
