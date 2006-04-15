# $Id$

SUBDIRS = test

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libquadrature.la

libquadrature_la_SOURCES = quadraturerules.cc
libquadrature_la_LIBADD = -lcommon

quadratureincludedir = $(includedir)/dune/quadrature
quadratureinclude_HEADERS = quadraturerules.hh

include $(top_srcdir)/am/global-rules
