# $Id: Makefile.am 2549 2006-04-28 13:09:28Z christi $

alugriddir = $(includedir)/dune/grid/alugrid/
alugrid_HEADERS = interfaces.hh defaultindexsets.hh

SUBDIRS = 2d 3d

noinst_LTLIBRARIES = libalu.la

# also contians items independent of ALUGrid
ALULIB = 3d/libalu3d.la

# construct only from sub-libraries
libalu_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libalu_la_SOURCES = dummy.cc
libalu_la_LIBADD = $(ALULIB)

include $(top_srcdir)/am/global-rules
