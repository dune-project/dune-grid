# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST= CMakeLists.txt dune.module

# don't follow the full GNU-standard
# we need automake 1.9 or newer
AUTOMAKE_OPTIONS = foreign 1.9

# All subdirectories that must (or can) be handled before the library is built
SUBDIRS =					\
        cmake                                   \
	doc					\
	dune					\
	m4
# The library subdirectory itself
SUBDIRS += lib
# All subdirectories that must (or can) be handled after the library is built
if EXTRA_UTILS
# Using automake conditionals ensures that these directories will always be in
# DIST_SUBDIRS, even if they are not always in SUBDIRS
SUBDIRS +=					\
	src
endif EXTRA_UTILS

# use configured compiler for "make distcheck"
# doxygen is difficult to convince to build in a special directory...
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" --with-dune-geometry="$(DUNE_GEOMETRY_ROOT)" CXX="$(CXX)" CC="$(CC)" --enable-parallel=@ENABLE_PARALLEL@ MPICC="$(MPICC)"

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/top-rules

# Distribute and install config.h.cmake
configdir = $(datadir)/dune-grid
dist_config_DATA = config.h.cmake
