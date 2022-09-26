# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# Module providing convenience methods for compile binaries with Alberta support.
#
# .. cmake_function:: add_dune_alberta_flags
#
#    .. cmake_param:: targets
#       :single:
#       :required:
#       :positional:
#
#       The targets to add the Alberta flags to.
#
#    .. cmake_param:: WORLDDIM
#       :single:
#
#       The dimension of the world space
#

# set HAVE_ALBERTA for config.h
set(HAVE_ALBERTA ${Alberta_FOUND})

# register all Alberta related flags
if(Alberta_FOUND)
  include(GridType)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM ASSERTION WORLDDIM == ALBERTA_DIM
    GRIDTYPE ALBERTAGRID DUNETYPE "Dune::AlbertaGrid< dimgrid >"
    HEADERS dune/grid/albertagrid.hh dune/grid/albertagrid/dgfparser.hh)
endif(Alberta_FOUND)

macro(add_dune_alberta_flags)
  if(Alberta_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(_ARG "" "WORLDDIM;GRIDDIM" "" ${ARGN})
    set(_ARG_TARGETS ${_ARG_UNPARSED_ARGUMENTS})

    # extract dimension arguments
    set(ALBERTA_DIM ${_ARG_WORLDDIM} ${_ARG_GRIDDIM} 0)
    list(SORT ALBERTA_DIM ORDER DESCENDING)
    list(GET ALBERTA_DIM 0 WORLDDIM)
    if(NOT WORLDDIM OR WORLDDIM EQUAL 0)
      message(WARNING "Alberta dimension not set. Please set it, e.g. use "
                      "add_dune_alberta_flags(WORLDDIM 2 <target>). "
                      "Falling back to dimension 2.")
      set(WORLDDIM 2)
    elseif(NOT ${WORLDDIM} IN_LIST ALBERTA_WORLD_DIMS)
      message(STATUS "Found ALBERTA_WORLD_DIMS=${ALBERTA_WORLD_DIMS}")
      message(FATAL_ERROR "There is no alberta library for dimension ${WORLDDIM}.")
    endif()

    # link to ALBERTA libraries
    foreach(_target ${_ARG_TARGETS})
      target_link_libraries(${_target} PUBLIC dunealbertagrid${WORLDDIM}d)
    endforeach(_target)
  endif(Alberta_FOUND)
endmacro(add_dune_alberta_flags)
