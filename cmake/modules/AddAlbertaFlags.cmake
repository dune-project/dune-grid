#
# Module providing convenience methods for compile binaries with Alberta support.
#
# Provides the following functions:
#
# add_dune_superlu_flags(target1 target2 ...)
#
# adds Alberta flags to the targets for compilation and linking
#
macro(add_dune_alberta_flags)
  if(ALBERTA_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_ALBERTA "OBJECT;SOURCE_ONLY;USE_GENERIC" "GRIDDIM;WORLDDIM" "" ${ARGN})

    if(ADD_ALBERTA_GRIDDIM AND NOT ADD_ALBERTA_WORLDDIM)
      set(ADD_ALBERTA_WORLDDIM ${ADD_ALBERTA_GRIDDIM})
    endif(ADD_ALBERTA_GRIDDIM AND NOT ADD_ALBERTA_WORLDDIM)

    if(ADD_ALBERTA_WORLDDIM AND NOT ADD_ALBERTA_GRIDDIM)
      set(ADD_ALBERTA_GRIDDIM ${ADD_ALBERTA_WORLDDIM})
    endif(ADD_ALBERTA_WORLDDIM AND NOT ADD_ALBERTA_GRIDDIM)

    if(NOT ADD_ALBERTA_WORLDDIM)
      message(WARNING "Alberta dimension not set. Please set it, e.g. use add_dune_alberta_flags(GRIDDIM 2 <target>). Falling back to dimension 2.")
      set(ADD_ALBERTA_WORLDDIM 2)
      set(ADD_ALBERTA_GRIDDIM ${ADD_ALBERTA_WORLDDIM})
    else(NOT ADD_ALBERTA_WORLDDIM)
      list(FIND ALBERTA_WORLD_DIMS ${ADD_ALBERTA_WORLDDIM} _index)
      if(_index EQUAL -1)
        message(FATAL_ERROR "There is no alberta library for dimension ${ADD_ALBERTA_GRIDDIM}.")
      endif(_index EQUAL -1)
    endif(NOT ADD_ALBERTA_WORLDDIM)
    if(ADD_ALBERTA_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${ALBERTA_INCLUDES})
    else()
      if(ADD_ALBERTA_OBJECT)
        # This a dune object libraries. Therefore we need to set the options
        # on the source files directly
        set(_all_sources "")
        foreach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
          get_property(_sources GLOBAL PROPERTY DUNE_LIB_${_target}_SOURCES)
          include_directories(${ALBERTA_INCLUDES})
          list(APPEND _all_sources ${_sources})
        endforeach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
        set(ADD_ALBERTA_UNPARSED_ARGUMENTS ${_all_sources}) #override unparsed arguments
      set(_prefix SOURCE)
      else(ADD_ALBERTA_OBJECT)
        # link to ALBERTA libraries
        foreach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} dunealbertagrid_${ADD_ALBERTA_GRIDDIM}d
            ${ALBERTA_${ADD_ALBERTA_GRIDDIM}D_LIB}
            dunegrid ${DUNE_LIBS} ${ALBERTA_UTIL_LIB} ${ALBERTA_EXTRA_LIBS})
        endforeach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
      set(_prefix TARGET)
      include_directories(${ALBERTA_INCLUDES})
      endif(ADD_ALBERTA_OBJECT)
    endif(ADD_ALBERTA_SOURCE_ONLY)

    replace_properties(${_prefix} ${ADD_ALBERTA_UNPARSED_ARGUMENTS}
      PROPERTY COMPILE_DEFINITIONS ENABLE_ALBERTA ENABLE_ALBERTA
      ALBERTA_DIM=.* ALBERTA_DIM=${ADD_ALBERTA_GRIDDIM}
      WORLDDIM=.* WORLDDIM=${ADD_ALBERTA_WORLDDIM})
    if(ADD_ALBERTA_USE_GENERIC)
      set_property(${_prefix} ${ADD_ALBERTA_UNPARSED_ARGUMENTS} APPEND
      PROPERTY COMPILE_DEFINITIONS  DUNE_ALBERTA_USE_GENERICGEOMETRY=1)
    endif(ADD_ALBERTA_USE_GENERIC)
  endif(ALBERTA_FOUND)
endmacro(add_dune_alberta_flags)
