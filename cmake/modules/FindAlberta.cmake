# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#[=======================================================================[.rst:
FindAlberta
-----------

Find Alberta, an Adaptive multiLevel finite element toolbox using Bisectioning
refinement and Error control by Residual Techniques for scientific Applications.
(see https://gitlab.mathematik.uni-stuttgart.de/ians-nmh/alberta/alberta3)

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` target:

``Alberta::AlbertaGrid_[n]d``
  Dimension dependent library

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``Alberta_FOUND``
  The Alberta library with all its dependencies is found

Cache Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``ENV{PKG_CONFIG_PATH}``
  Directory containing the `alberta-grid_[n]d.pc` pkg-config files.
  An environmental variable to influence the search procedure of pkg-config
  for finding Alberta.

``Alberta_FIND_QUIETLY``
  If set to `ON` do not print detailed information during search for
  Alberta using pkg-config. This variable is automatically if `find_package`
  if invoked with flag `QUIET`.

``ALBERTA_MAX_WORLD_DIM``
  Maximal world dimension to check for Alberta library. Default: 3.

#]=======================================================================]

# text for feature summary
include(FeatureSummary)
set_package_properties("Alberta" PROPERTIES
  DESCRIPTION "An adaptive hierarchical finite element toolbox and grid manager")

set(ALBERTA_MAX_WORLD_DIM "3" CACHE STRING "Maximal world dimension to check for Alberta library.")

set(ALBERTA_WORLD_DIMS)
set(ALBERTA_GRID_VERSION)
set(ALBERTA_GRID_PREFIX)

# search for Alberta using pkg-config
find_package(PkgConfig)
if(PkgConfig_FOUND)
  foreach(dim RANGE 1 ${ALBERTA_MAX_WORLD_DIM})
    if(Alberta_FIND_VERSION)
      set(ALBERTA_GRID_DIM_MODULE "alberta-grid_${dim}d>=${Alberta_FIND_VERSION}")
    else()
      set(ALBERTA_GRID_DIM_MODULE "alberta-grid_${dim}d")
    endif()

    if(Alberta_FIND_QUIETLY)
      pkg_check_modules(Alberta${dim}d ${ALBERTA_GRID_DIM_MODULE} QUIET IMPORTED_TARGET GLOBAL)
    else()
      pkg_check_modules(Alberta${dim}d ${ALBERTA_GRID_DIM_MODULE} IMPORTED_TARGET GLOBAL)
    endif()
    if(Alberta${dim}d_FOUND)
      list(APPEND ALBERTA_WORLD_DIMS ${dim})
      set(ALBERTA_GRID_VERSION ${Alberta${dim}d_VERSION})
      set(ALBERTA_GRID_PREFIX ${Alberta${dim}d_PREFIX})
    endif()
  endforeach(dim)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("Alberta"
  REQUIRED_VARS
    ALBERTA_GRID_PREFIX PkgConfig_FOUND
  VERSION_VAR
    ALBERTA_GRID_VERSION
  FAIL_MESSAGE "Could NOT find Alberta (set PKG_CONFIG_PATH to include the location of the alberta-grid_[n]d.pc files)"
)

if(Alberta_FOUND)
  foreach(dim ${ALBERTA_WORLD_DIMS})
    if(NOT Alberta::AlbertaGrid_${dim}d)
      add_library(Alberta::AlbertaGrid_${dim}d ALIAS PkgConfig::Alberta${dim}d)
    endif()
  endforeach(dim)
endif()
