# If dune-uggrid was found this module adds the dgf magic to config.h
# and makes add_dune_ug_flags available.
#
# .. cmake_function:: add_dune_ug_flags
#
#    .. cmake_param:: targets
#       :single:
#       :positional:
#       :required:
#
#       The targets to add the UG flags to.
#
#    .. cmake_param:: SOURCE_ONLY
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are source files.
#
#    .. cmake_param:: OBJECT
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are object libraries.
#

# Add dgf magic to config.h and register flags
if(dune-uggrid_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "${UG_DEFINITIONS}"
                              INCLUDE_DIRS "${UG_INCLUDES}"
                              LIBRARIES "dunegrid;${UG_LIBRARIES};${DUNE_LIBS}")

  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID
    ASSERTION "GRIDDIM == WORLDDIM"
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS "dune/grid/uggrid.hh" "dune/grid/io/file/dgfparser/dgfug.hh")
endif()

# Add flags to targets
function(add_dune_ug_flags)
  if(dune-uggrid_FOUND)
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT;NO_LINK_DUNEGRID" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
    else()
      if(NOT ADD_UG_OBJECT)
        foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
          if(NOT ADD_UG_NO_LINK_DUNEGRID)
            target_link_libraries(${_target} PUBLIC dunegrid)
          endif()
          target_link_libraries(${_target}
            PUBLIC ${UG_LIBRARIES} ${DUNE_LIBS})
        endforeach()
      endif()
      set(_prefix TARGET)
    endif()

    if(ADD_UG_OBJECT)
      set(_object OBJECT)
    else()
      set(_object)
    endif()

    # Add compiler flags
    include_directories(${UG_INCLUDES})
    set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
      APPEND PROPERTY
      COMPILE_DEFINITIONS ${UG_DEFINITIONS})
    if(UG_PARALLEL)
      # Add mpi flags.
      add_dune_mpi_flags(${ADD_UG_UNPARSED_ARGUMENTS} ${_source_only} ${_object})
    endif()
  endif()
endfunction(add_dune_ug_flags)
