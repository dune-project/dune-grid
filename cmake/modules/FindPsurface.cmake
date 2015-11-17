# .. cmake_module::
#
#    Module that checks whether psurface is available and usable.
#
#    Variables used by this module which you may want to set:
#
#    :ref:`PSURFACE_ROOT`
#       Path list to search for psurface
#
#    Sets the follwing variables:
#
#    :code:`PSURFACE_FOUND`
#       True if PSurface available and usable.
#
#    :code:`PSURFACE_INCLUDE_DIRS`
#       Path to the PSurface include dirs.
#
#    :code:`PSURFACE_LIBRARIES`
#       Name to the PSurface library.
#
# .. cmake_variable:: PSURFACE_ROOT
#
#    You may set this variable to have :ref:`FindPsurface` look
#    for the Psurface package in the given path before inspecting
#    system paths.
#

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

set(PSURFACE_WITH_VERSION "psurface 1.3" CACHE STRING
  "Human readable string containing psurface version information.")

# Try to find psurface with pkg-config (for psurface 2.0 or newer)
include(FindPkgConfig)
pkg_search_module(PKG_PSURFACE psurface)

if(NOT PKG_PSURFACE_FOUND)
  # first only at positions given by the user
  find_file(PATH_PKG_PSURFACE
    NAMES "psurface.pc"
    PATHS ${PSURFACE_ROOT}
    PATH_SUFFIXES lib/pkgconfig lib32/pkgconfig lib64/pkgconfig
    NO_DEFAULT_PATH)
  # including default paths
  find_file(PATH_PKG_PSURFACE
    NAMES "psurface.pc"
    PATH_SUFFIXES lib/pkgconfig lib32/pkgconfig lib64/pkgconfig)

  # try again with path temporarilly added to PKG_CONFIG_PATH
  set(REM_PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
  get_filename_component(DIR_PKG_PSURFACE "${PATH_PKG_PSURFACE}" PATH)
  set(ENV{PKG_CONFIG_PATH} "${PSURFACE_ROOT}:${DIR_PKG_PSURFACE}:${PKG_CONFIG_PATH}")
  pkg_check_modules(PKG_PSURFACE psurface)
  set(ENV{PKG_CONFIG_PATH} REM_PKG_CONFIG_PATH)
endif()

if(PKG_PSURFACE_FOUND)
  set(HAVE_PSURFACE_2_0 1)
  set(PSURFACE_WITH_VERSION "psurface >= 2.0" CACHE STRING
    "Human readable string containing psurface version information.")
endif()
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
    "Determining location of ${PSURFACE_WITH_VERSION} succeeded:\n"
    "Include directory: ${PSURFACE_INCLUDE_DIRS}\n"
    "Library directory: ${PSURFACE_LIBRARIES}\n\n")
  set(PSURFACE_DUNE_COMPILE_FLAGS "-I${PSURFACE_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling psurface programs")
  set(PSURFACE_DUNE_LIBRARIES ${PSURFACE_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking psurface programs")
else(PSURFACE_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of psurface failed:\n"
    "Include directory: ${PSURFACE_INCLUDE_DIRS}\n"
    "Library directory: ${PSURFACE_LIBRARIES}\n\n")
endif(PSURFACE_FOUND)

# set HAVE_PSURFACE for config.h
set(HAVE_PSURFACE ${PSURFACE_FOUND})

# register all psurface related flags
if(PSURFACE_FOUND)
  dune_register_package_flags(INCLUDE_DIRS "${PSURFACE_INCLUDE_DIRS}"
                              LIBRARIES "${PSURFACE_LIBRARIES}")
endif()
