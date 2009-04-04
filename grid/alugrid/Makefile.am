# $Id: Makefile.am 2549 2006-04-28 13:09:28Z christi $

alugriddir = $(includedir)/dune/grid/alugrid/
alugrid_HEADERS = dynamiccodimsubindexid.hh

SUBDIRS = 2d 3d

noinst_LTLIBRARIES = libalu.la

# also contians items independent of ALUGrid
ALULIB = 3d/libalu3d.la

# construct only from sub-libraries
libalu_la_SOURCES =
libalu_la_LIBADD = $(ALULIB)

include $(top_srcdir)/am/global-rules
