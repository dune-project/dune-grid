#
# Module providing convenience methods for compile binaries with Grape support.
#
# Provides the following functions:
#
# add_dune_superlu_flags(target1 target2 ...)
#
# adds Grape flags to the targets for compilation and linking
#
function(add_dune_grape_flags)
  if(GRAPE_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_GRAPE "OBJECT;SOURCE_ONLY" "" "" ${ARGN})
    if(ADD_GRAPE_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${GRAPE_INCLUDE_DIRS})
    else(ADD_GRAPE_SOURCE_ONLY)
      if(NOT ADD_GRAPE_OBJECT)
        foreach(_target ${ADD_GRAPE_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${GRAPE_LIBARIES})
        endforeach(_target ${ADD_GRAPE_UNPARSED_ARGUMENTS})
      endif(NOT ADD_GRAPE_OBJECT)
      set(_prefix TARGET)

      set_property(${_prefix}  ${ADD_GRAPE_UNPARSED_ARGUMENTS} APPEND
        PROPERTY
        COMPILE_DEFINITIONS -DENABLE_GRAPE=1)
      include_directories(${GRAPE_INCLUDE_DIRS})
    endif(ADD_GRAPE_SOURCE_ONLY)
  endif(GRAPE_FOUND)
endfunction(add_dune_grape_flags)
