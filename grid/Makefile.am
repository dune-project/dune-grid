# $Id$

SUBDIRS = albertagrid alugrid common onedgrid sgrid uggrid utility yaspgrid \
	test io
noinst_LTLIBRARIES = libgrid.la

# always part of lib:
ONEDLIB = onedgrid/libonedgrid.la
COMMONLIB = common/libgridcommon.la io/libio.la

# conditional parts:
if UG
  UGLIB = uggrid/libuggrid.la
endif

if ALUGRID
# also contains items independent of ALUGrid
ALULIB = alugrid/libalu.la
endif

if ALBERTA
ALBERTALIB = albertagrid/libalberta.la
endif

# construct only from sub-libraries
libgrid_la_SOURCES =
libgrid_la_LIBADD = $(ONEDLIB) $(UGLIB) $(ALULIB) $(ALBERTALIB) \
  $(COMMONLIB) $(DUNE_COMMON_LIBS)

# These libraries do not need other libraries besides lcommon
griddir = $(includedir)/dune/grid
grid_HEADERS = albertagrid.hh \
               alugrid.hh \
               sgrid.hh \
               uggrid.hh \
               yaspgrid.hh \
               onedgrid.hh

include $(top_srcdir)/am/global-rules
