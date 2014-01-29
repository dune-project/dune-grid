#
# Module that checks whether psurface is available and usable.
#
# Variables used by this module which you may want to set:
# PSURFACE_ROOT         Path list to search for psurface
#
# Sets the follwing variable:
#
# PSURFACE_FOUND          True if PSurface available and usable.
# PSURFACE_INCLUDE_DIRS   Path to the PSurface include dirs.
# PSURFACE_LIBRARIES      Name to the PSurface library.
#
message(AUTHOR_WARNING "TODO: Implement Amiramesh support for PSurface test")

# look for header files, only at positions given by the user
find_path(PSURFACE_INCLUDE_DIR
  NAMES "psurface/PSurface.h"
  PATHS ${PSURFACE_PREFIX} ${PSURFACE_ROOT}
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH
)

# look for header files, including default paths
find_path(PSURFACE_INCLUDE_DIR
  NAMES "psurface/PSurface.h"
  PATH_SUFFIXES "include"
)

# look for library, only at positions given by the user
find_library(PSURFACE_LIBRARY
  NAMES "psurface"
  PATHS ${PSURFACE_PREFIX} ${PSURFACE_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(PSURFACE_LIBRARY
  NAMES "psurface"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)

# check version specific macros
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND
# if the searches above were not successful
# Without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
# Please set them or make sure they are set and tested correctly in the CMake files:"
#
if(PSURFACE_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${PSURFACE_INCLUDE_DIR})
endif(PSURFACE_INCLUDE_DIR)
if(PSURFACE_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${PSURFACE_LIBRARY})
endif(PSURFACE_LIBRARY)

# Try to link to the library (for libpsurface-1.3 and newer)
CHECK_CXX_SOURCE_COMPILES("
#include <psurface/PSurface.h>
int main(void)
{
  psurface::PSurface<2,double> foo;
  return 0;
}"
PSURFACE_MIN_VERSION_1_3)

cmake_pop_check_state()

if(PSURFACE_MIN_VERSION_1_3)
  set(PSURFACE_WITH_VERSION "psurface >= 1.3" CACHE STRING
    "Human readable string containing psurface version information.")
  set(PSURFACE_NAMESPACE "psurface::")
else()
  set(PSURFACE_WITH_VERSION "psurface <= 1.2" CACHE STRING
    "Human readable string containing psurface version information.")
  set(PSURFACE_NAMESPACE "")
endif(PSURFACE_MIN_VERSION_1_3)

# Try to find psurface with pkg-config (for psurface 2.0 or newer)
include(FindPkgConfig)
pkg_search_module(PKG_PSURFACE psurface)
if(${PKG_PSURFACE_FOUND})
  set(HAVE_PSURFACE_2_0 1)
  set(PSURFACE_WITH_VERSION "psurface >= 2.0" CACHE STRING
    "Human readable string containing psurface version information.")
endif(${PKG_PSURFACE_FOUND})
# re-set PKG_CONFIG_PATH
set(PKG_CONFIG_PATH ${PKG_CONFIG_PATH_STORE})

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "psurface"
  DEFAULT_MSG
  PSURFACE_INCLUDE_DIR
  PSURFACE_LIBRARY
)

mark_as_advanced(PSURFACE_INCLUDE_DIR PSURFACE_LIBRARY PKG_PSURFACE_FOUND)

# if both headers and library are found, store results
if(PSURFACE_FOUND)
  set(PSURFACE_INCLUDE_DIRS ${PSURFACE_INCLUDE_DIR})
  set(PSURFACE_LIBRARIES    ${PSURFACE_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ${PSURFACE_WITH_VERSION} succeded:\n"
    "Include directory: ${PSURFACE_INCLUDE_DIRS}\n"
    "Library directory: ${PSURFACE_LIBRARIES}\n\n")
  set(PSURFACE_DUNE_COMPILE_FLAGS "-I${PSURFACE_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling psurface programs")
  set(PSURFACE_DUNE_LIBRARIES ${PSURFACE_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking psurface programs")
else(PSURFACE_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of psurface failed:\n"
    "Include directory: ${PSURFACE_INCLUDE_DIRS}\n"
    "Library directory: ${PSURFACE_LIBRARIES}\n\n")
endif(PSURFACE_FOUND)

# set HAVE_PSURFACE for config.h
set(HAVE_PSURFACE ${PSURFACE_FOUND})

#add all psurface related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_psurface_flags
if(PSURFACE_FOUND)
  foreach(dir ${PSURFACE_INCLUDE_DIRS})
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()
