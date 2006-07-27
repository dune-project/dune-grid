# $Id$

# don't follow the full GNU-standard
# need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = grid doc lib

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dune-grid.pc

# use configured compiler for "make distcheck"
# doxygen is difficult to convince to build in a special directory...
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" CXX="$(CXX)" CC="$(CC)" DOXYGEN="true"

include $(top_srcdir)/am/global-rules
