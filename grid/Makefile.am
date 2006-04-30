# $Id$

SUBDIRS = albertagrid alugrid common onedgrid sgrid uggrid utility yaspgrid \
	test

noinst_LTLIBRARIES = libgrid.la

# always part of lib:
ONEDLIB = onedgrid/libonedgrid.la
COMMONLIB = common/libgridcommon.la io/libio.la

# conditional parts:
if UG
  UGLIB = uggrid/libuggrid.la
endif

# also contians items independent of ALUGrid
ALULIB = alugrid/libalu.la

# construct only from sub-libraries
libgrid_la_SOURCES =
libgrid_la_LIBADD = $(ONEDLIB) $(UGLIB) $(ALULIB) $(COMMONLIB) -lcommon

griddir = $(includedir)/dune/grid
grid_HEADERS = albertagrid.hh \
               alu3dgrid.hh \
               alugrid.hh \
               sgrid.hh \
               uggrid.hh \
               yaspgrid.hh \
               onedgrid.hh

include $(top_srcdir)/am/global-rules
