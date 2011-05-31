# $Id$

SUBDIRS = genericgeometry common alugrid onedgrid sgrid uggrid utility \
          yaspgrid io albertagrid geometrygrid test
noinst_LTLIBRARIES = libgrid.la

# always part of lib:
ONEDLIB = onedgrid/libonedgrid.la
COMMONLIB = common/libgridcommon.la
DGFPARSERLIB = io/file/dgfparser/libdgfparser.la

# conditional parts:
if UG
  UGLIB = uggrid/libuggrid.la
endif

if AMIRAMESH
if PSURFACE
  AMIRAMESHLIB = io/file/amiramesh/libamiramesh.la
endif
endif

# the ALUGrid lib contains items that are needed
# in any case, so no dependency from ALUGRID flag
ALULIB = alugrid/libalu.la

# construct only from sub-libraries
libgrid_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libgrid_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libgrid_la_LIBADD = $(ONEDLIB) $(UGLIB) $(ALULIB) \
	$(COMMONLIB)				\
	$(DGFPARSERLIB)				\
        $(AMIRAMESHLIB)                         \
	$(DUNE_LIBS)

# These libraries do not need other libraries besides lcommon
griddir = $(includedir)/dune/grid
grid_HEADERS = albertagrid.hh \
               alugrid.hh \
               sgrid.hh \
               uggrid.hh \
               yaspgrid.hh \
               onedgrid.hh \
	       geometrygrid.hh

include $(top_srcdir)/am/global-rules
