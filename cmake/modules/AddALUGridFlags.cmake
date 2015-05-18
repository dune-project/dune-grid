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
    else()
      if(ADD_ALU_OBJECT)
        set(_object OBJECT)
      else(ADD_ALU_OBJECT)
        foreach(_target ${ADD_ALU_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target}
            dunegrid ${ALUGRID_LIBRARIES} ${METIS_LIBRARIES} ${DUNE_LIBS})
        endforeach(_target ${ADD_ALU_UNPARSED_ARGUMENTS})
      endif(ADD_ALU_OBJECT)
      set(_prefix TARGET)
    endif()

    include_directories(${ALUGRID_INCLUDES})
    set_property(${_prefix} ${ADD_ALU_UNPARSED_ARGUMENTS}
      APPEND PROPERTY
      COMPILE_DEFINITIONS ENABLE_ALUGRID)
    # add linker arguments
    if(NOT (ADD_ALU_SOURCE_ONLY OR ADD_ALU_OBJECT))
      set_property(${_prefix} ${ADD_ALU_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        LINK_LIBRARIES dunegrid ${ALUGRID_LIBRARIES} ${METIS_LIBRARIES} ${DUNE_LIBS})
    endif(NOT (ADD_ALU_SOURCE_ONLY OR ADD_ALU_OBJECT))
    if(HAVE_ALUGRID_METIS AND (NOT ADD_ALU_SOURCE_ONLY))
      if(PARMETIS_FOUND)
        include(AddParMETISFlags)
        add_dune_parmetis_flags(${ADD_ALU_UNPARSED_ARGUMENTS})
      else()
        include(AddMETISFlags)
        add_dune_metis_flags(${ADD_ALU_UNPARSED_ARGUMENTS})
      endif()
      add_dune_mpi_flags(${ADD_ALU_UNPARSED_ARGUMENTS} ${_source_only} ${_object})
    endif(HAVE_ALUGRID_METIS AND (NOT ADD_ALU_SOURCE_ONLY))
  endif(ALUGRID_FOUND)
endfunction(add_dune_alugrid_flags)
