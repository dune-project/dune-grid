#[=======================================================================[.rst:
FindAlberta
-----------

Find Alberta, an Adaptive multiLevel finite element toolbox using Bisectioning
refinement and Error control by Residual Techniques for scientific Applications.
(see https://gitlab.mathematik.uni-stuttgart.de/ians-nmh/alberta/alberta3)

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` target:

``Alberta::Alberta``
  The libraries, flags, and includes to use for Alberta, if found.

``Alberta::Alberta_XD``
  Dimension dependent library

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``Alberta_FOUND``
  The Alberta library with all its dependencies is found

Cache Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``Alberta_ROOT``
  Base directory of Alberta installation

``ALBERTA_EXTRA_LIBS``
  Libraries needed to link Alberta, e.g. m

#]=======================================================================]

# text for feature summary
include(FeatureSummary)
set_package_properties("Alberta" PROPERTIES
  DESCRIPTION "An adaptive hierarchical finite element toolbox and grid manager")

set(ALBERTA_EXTRA_LIBS "m" CACHE FILEPATH "Extra libraries needed by alberta for linking.")

# look for header alberta/alberta.h
find_path(ALBERTA_INCLUDE_DIR "alberta/alberta.h")

# check version, Alberta 2 or 3
set(ALBERTA_VERSION 2.0)
set(DUNE_ALBERTA_VERSION 0x200)
find_file(ALBERTA_HEADER "alberta/alberta.h"
  HINTS ${ALBERTA_INCLUDE_DIR}
  NO_DEFAULT_PATH)
if(ALBERTA_HEADER)
  cmake_push_check_state()
  list(APPEND CMAKE_REQUIRED_DEFINITIONS -DDIM_OF_WORLD=3 -DDEL_INDEX=0)
  list(APPEND CMAKE_REQUIRED_INCLUDES ${ALBERTA_INCLUDE_DIR})

  include(CheckStructHasMember)
  check_struct_has_member("struct el_info" wall_bound ${ALBERTA_HEADER} ALBERTA_IS_VERSION_3)
  if(ALBERTA_IS_VERSION_3)
    set(ALBERTA_VERSION 3.0)
    set(DUNE_ALBERTA_VERSION 0x300)
  endif()
  cmake_pop_check_state()
endif()
unset(ALBERTA_HEADER CACHE)

# check for ltdl
find_library(ALBERTA_LTDL_LIB "ltdl")
if(ALBERTA_LTDL_LIB)
  list(APPEND ALBERTA_EXTRA_LIBS "ltdl")
endif()

# look for alberta utilitiy libraries
find_library(ALBERTA_UTIL_LIB
  NAMES alberta_util alberta_utilities)


# check for which dimensions are supported Alberta installation
cmake_push_check_state()
list(APPEND CMAKE_REQUIRED_LIBRARIES "${ALBERTA_EXTRA_LIBS};${ALBERTA_UTIL_LIB}")

foreach(dim RANGE 1 9)
  find_library(ALBERTA_${dim}D_LIB alberta_${dim}d)
  if(ALBERTA_${dim}D_LIB)
    check_library_exists(${ALBERTA_${dim}D_LIB} mesh_traverse "" ALBERTA_${dim}D_LIB_FOUND)
    if(ALBERTA_${dim}D_LIB_FOUND)
      list(APPEND ALBERTA_WORLD_DIMS ${dim})
    endif()
  endif()
  mark_as_advanced(ALBERTA_${dim}D_LIB)
endforeach(dim)

cmake_pop_check_state()
message(STATUS "Found alberta libraries for dimensions ${ALBERTA_WORLD_DIMS}")

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("Alberta"
  REQUIRED_VARS
    ALBERTA_INCLUDE_DIR ALBERTA_UTIL_LIB ALBERTA_WORLD_DIMS
  VERSION_VAR
    ALBERTA_VERSION
)

mark_as_advanced(ALBERTA_INCLUDE_DIR ALBERTA_UTIL_LIB ALBERTA_WORLD_DIMS)

if(Alberta_FOUND AND NOT TARGET Alberta::Alberta)
  add_library(Alberta::Alberta UNKNOWN IMPORTED)
  set_target_properties(Alberta::Alberta PROPERTIES
    IMPORTED_LOCATION ${ALBERTA_UTIL_LIB}
    INTERFACE_INCLUDE_DIRECTORIES ${ALBERTA_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES "${ALBERTA_EXTRA_LIBS}")

  foreach(dim RANGE 1 9)
    if(ALBERTA_${dim}D_LIB AND NOT Alberta::Alberta_${dim}D)
      add_library(Alberta::Alberta_${dim}D UNKNOWN IMPORTED)
      set_target_properties(Alberta::Alberta_${dim}D PROPERTIES
        IMPORTED_LOCATION ${ALBERTA_${dim}D_LIB}
        INTERFACE_COMPILE_DEFINITIONS ALBERTA_DIM=${dim})
    endif()
  endforeach(dim)
endif()
