# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard
# need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = dune doc lib m4

# use configured compiler for "make distcheck"
# doxygen is difficult to convince to build in a special directory...
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" CXX="$(CXX)" CC="$(CC)" --enable-parallel=@ENABLE_PARALLEL@ MPICC="$(MPICC)"

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/top-rules
