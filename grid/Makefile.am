# $Id$

SUBDIRS = albertagrid alugrid common onedgrid sgrid uggrid utility yaspgrid \
	test io
noinst_LTLIBRARIES = libgrid.la

# always part of lib:
ONEDLIB = onedgrid/libonedgrid.la
COMMONLIB = common/libgridcommon.la

# conditional parts:
if UG
  UGLIB = uggrid/libuggrid.la
endif

if AMIRAMESH
if UG
  AMIRALIB = io/file/amiramesh/libugamiramesh.la
endif
endif

# the ALUGrid lib contains items that are needed
# in any case, so no dependency from ALUGRID flag
ALULIB = alugrid/libalu.la

# construct only from sub-libraries
libgrid_la_SOURCES =
libgrid_la_LIBADD = $(ONEDLIB) $(UGLIB) $(ALULIB) \
  $(COMMONLIB) $(DUNE_COMMON_LIBS) $(AMIRALIB)

# These libraries do not need other libraries besides lcommon
griddir = $(includedir)/dune/grid
grid_HEADERS = albertagrid.hh \
               alugrid.hh \
               sgrid.hh \
               uggrid.hh \
               yaspgrid.hh \
               onedgrid.hh

include $(top_srcdir)/am/global-rules
