#[=======================================================================[.rst:
FindPsurface
------------

Find the Psurface library, a C++ library that handles piecewise linear
bijections between triangulated surfaces.
(see https://github.com/psurface/psurface)

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` target:

``Psurface::Psurface``
  The libraries, flags, and includes to use for Psurface, if found.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``Psurface_FOUND``
  The Psurface library with all its dependencies is found

Cache Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``PSURFACE_INCLUDE_DIR``
  Include directory where the PSurface.h is found.

``PSURFACE_LIBRARY``
  Full path to the psurface library

Environmental Variables
^^^^^^^^^^^^^^^^^^^^^^^

The find module uses PkgConfig in case the header and library are not found
directly. By adding the directory where to find the `psurface.pc` file,
PkgConfig can be directed to the library:

``ENV{PKG_CONFIG_PATH}``
  Environmental variable of search-paths for PkgConfig.

#]=======================================================================]

# text for feature summary
include(FeatureSummary)
set_package_properties("Psurface" PROPERTIES
  DESCRIPTION "Piecewise linear bijections between triangulated surfaces"
)

# look for header files
find_path(PSURFACE_INCLUDE_DIR "PSurface.h"
  PATH_SUFFIXES psurface)

# look for the psurface library
find_library(PSURFACE_LIBRARY "psurface")

if(NOT PSURFACE_INCLUDE_DIR OR NOT PSURFACE_LIBRARY)
  # Try to find psurface with pkg-config (for psurface 2.0 or newer)
  include(FindPkgConfig)
  pkg_search_module(PKG_PSURFACE psurface IMPORTED_TARGET)

  if(NOT PKG_PSURFACE_FOUND)
    find_file(PATH_PKG_PSURFACE psurface.pc
      PATH_SUFFIXES lib/pkgconfig lib32/pkgconfig lib64/pkgconfig)

    # try again with path temporarilly added to PKG_CONFIG_PATH
    set(REM_PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
    get_filename_component(DIR_PKG_PSURFACE "${PATH_PKG_PSURFACE}" PATH)
    set(ENV{PKG_CONFIG_PATH} "${DIR_PKG_PSURFACE}:${PKG_CONFIG_PATH}")
    pkg_check_modules(PKG_PSURFACE psurface IMPORTED_TARGET)
    set(ENV{PKG_CONFIG_PATH} REM_PKG_CONFIG_PATH)
  endif()

  if(PKG_PSURFACE_FOUND)
    set(PSURFACE_INCLUDE_DIR ${PKG_PSURFACE_INCLUDE_DIRS})
    set(PSURFACE_LIBRARY ${PKG_PSURFACE_LINK_LIBRARIES})
  endif()
endif()

mark_as_advanced(PSURFACE_INCLUDE_DIR PSURFACE_LIBRARY)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("Psurface"
  REQUIRED_VARS
    PSURFACE_INCLUDE_DIR
    PSURFACE_LIBRARY
)

if(Psurface_FOUND AND NOT TARGET Psurface::Psurface)
  if(PKG_PSURFACE_FOUND)
    add_library(Psurface::Psurface ALIAS PkgConfig::psurface)
  else()
    add_library(Psurface::Psurface UNKNOWN IMPORTED)
    set_target_properties(Psurface::Psurface PROPERTIES
      IMPORTED_LOCATION ${PSURFACE_LIBRARY}
      INTERFACE_INCLUDE_DIRECTORIES ${PSURFACE_INCLUDE_DIR}
    )
  endif()
endif()
