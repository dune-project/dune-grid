# This module first provides macros for setting the necessary flags
# when compiling with dune-uggrid support.

# Add flags to targets
function(add_dune_ug_flags)
  if(UG_FOUND)
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT;NO_LINK_DUNEGRID" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
    else()
      if(NOT ADD_UG_OBJECT)
        foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
          if(NOT ADD_UG_NO_LINK_DUNEGRID)
            target_link_libraries(${_target} dunegrid)
          endif()
          target_link_libraries(${_target}
            ${UG_LIBRARIES} ${DUNE_LIBS})
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
      COMPILE_DEFINITIONS ENABLE_UG)
    if (UG_VERSION VERSION_GREATER 3.13.0
        OR UG_VERSION VERSION_EQUAL 3.13.0)
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        COMPILE_DEFINITIONS UG_USE_NEW_DIMENSION_DEFINES)
    endif()
    # Add linker arguments
    if(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
      if(NOT ADD_UG_NO_LINK_DUNEGRID)
        set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
          APPEND PROPERTY
          LINK_LIBRARIES dunegrid)
      endif()
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        LINK_LIBRARIES ${UG_LIBRARIES} ${DUNE_LIBS})
    endif(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
    if(UG_PARALLEL STREQUAL "yes")
      # Add modelp
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        COMPILE_DEFINITIONS ModelP)
      # Add mpi flags.
      add_dune_mpi_flags(${ADD_UG_UNPARSED_ARGUMENTS} ${_source_only} ${_object})
    endif()
  endif()
endfunction(add_dune_ug_flags)
