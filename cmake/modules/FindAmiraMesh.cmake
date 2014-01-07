#
# Module that checks whether AmiraMesh is available and usable.
#
# Variables used by this module which you may want to set:
# AMIRAMESH_ROOT            Path list to search for AmiraMesh
#
# Sets the follwing variable:
#
# AMIRAMESH_FOUND           True if AmiraMesh available and usable.
# AMIRAMESH_INCLUDE_DIRS    Path to the AmiraMesh include dirs.
# AMIRAMESH_LIBRARIES       Name to the AmiraMesh library.
#

find_path(AMIRAMESH_INCLUDE_DIR
  NAMES "amiramesh/AmiraMesh.h"
  PATHS ${AMIRAMESH_PREFIX} ${AMIRAMESH_ROOT}
  PATH_SUFFIXES "include"
)

find_library(AMIRAMESH_LIBRARY
  NAMES libamiramesh.a amiramesh
  PATHS ${AMIRAMESH_PREFIX} ${AMIRAMESH_ROOT}
  PATH_SUFFIXES "lib"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "AmiraMesh"
  DEFAULT_MSG
  AMIRAMESH_INCLUDE_DIR
  AMIRAMESH_LIBRARY
)

mark_as_advanced(AMIRAMESH_INCLUDE_DIR AMIRAMESH_LIBRARY AMIRAMESH_FOUND)

set(HAVE_AMIRAMESH ${AMIRAMESH_FOUND})

if(AMIRAMESH_FOUND)
  set(AMIRAMESH_INCLUDE_DIRS "${AMIRAMESH_INCLUDE_DIR}")
  set(AMIRAMESH_LIBRARIES "${AMIRAMESH_LIBRARY}")
  set(AMIRAMESH_COMPILE_FLAGS "-I${AMIRAMESH_INCLUDE_DIR}")
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${AMIRAMESH_INCLUDE_DIR}")
endif()
