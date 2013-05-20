#
# Module providing convenience methods for compile binaries with ALUGrid support.
#
# Provides the following functions:
#
# add_dune_alugrid_flags(target1 target2 ...)
#
# adds ALUGrid flags to the targets for compilation and linking
#
function(add_dune_alugrid_flags )
 if(ALUGRID_FOUND)
    cmake_parse_arguments(ADD_ALU "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_ALU_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${ALUGRID_INCLUDES})
    else()
      if(ADD_ALU_OBJECT)
        set(_object OBJECT)
      else(ADD_ALU_OBJECT)
        foreach(_target ${ADD_ALU_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} dunegrid ${ALUGRID_LIBRARIES} ${METIS_LIBRARIES} dunegeometry dunecommon)
        endforeach(_target ${ADD_ALU_UNPARSED_ARGUMENTS})
      endif(ADD_ALU_OBJECT)
      set(_prefix TARGET)
      include_directories(${ALUGRID_INCLUDES})
    endif()

    set_property(${_prefix} ${ADD_ALU_UNPARSED_ARGUMENTS}
      APPEND PROPERTY COMPILE_DEFINITIONS ENABLE_ALUGRID)
    if(NOT (ADD_ALU_SOURCE_ONLY OR ADD_ALU_OBJECT))
      set_property(${_prefix} ${ADD_ALU_UNPARSED_ARGUMENTS}
        APPEND PROPERTY LINK_LIBRARIES  dunegrid ${ALUGRID_LIBRARIES}  ${METIS_LIBRARIES} dunegeometry dunecommon)
    endif(NOT (ADD_ALU_SOURCE_ONLY OR ADD_ALU_OBJECT))
    if(HAVE_ALUGRID_PARALLEL_H)
      add_dune_mpi_flags(${ADD_ALU_UNPARSED_ARGUMENTS} ${_source_only} ${_object})
    endif(HAVE_ALUGRID_PARALLEL_H)
  endif(ALUGRID_FOUND)
endfunction(add_dune_alugrid_flags)
