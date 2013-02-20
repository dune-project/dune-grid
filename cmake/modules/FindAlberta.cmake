macro(_dune_set_alberta val)
  set(ALBERTA_FOUND ${val})
  set(HAVE_ALBERTA $val})
endmacro(_dune_set_alberta val)

set(ALBERTA_LIBCHECK ON CACHE BOOL "Whether to try to link against libalberta_Nd.")
set(ALBERTA_DIR "" CACHE FILEPATH "Root directory of Alberta installation.")
set(ALBERTA_EXTRA_LIBS "" CACHE FILEPATH "Extra libraries needed by alberta for linking.")

find_path(ALBERTA_INCLUDE_DIR name alberta/alberta.h PATHS ${ALBERTA_DIR}
  PATH_SUFFIXES alberta NO_DEFAULT_PATH DOC "Include path of Alberta")
find_path(ALBERTA_INCLUDE_DIR name alberta/alberta.h PATHS /usr/local /opt
  PATH_SUFFIXES alberta)

set(ALBERTA_INCLUDES ${ALBERTA_INCLUDE_DIR})
set(ALBERTA_VERSION 2.0)

cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} -DDIM_OF_WORLD=3 -DDEL_INDEX=0)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALBERTA_INCLUDE_DIR})
check_include_files(alberta/alberta.h ALBERTA_FOUND)
if(ALBERTA_FOUND)
  include(CheckStructHasMember)
  check_struct_has_member ("struct el_info" wall_bound alberta/alberta.h ALBERTA_IS_VERSION_3)
  if(!ALBERTA_IS_VERSION_3)
    message(WARNING "Alberta version 3 not found. Checking for version 2.")
  endif(!ALBERTA_IS_VERSION_3)
else(ALBERTA_FOUND)
  message(FATAL_ERROR "Error checking alberta.h")
endif(ALBERTA_FOUND)

find_library(ALBERTA_UTIL_LIB alberta_util PATH ${ALBERTA_DIR} PATH_SUFFIXES lib lib32 lib64
  NO_DEFAULT_PATH)
find_library(ALBERTA_UTIL_LIB alberta_util PATH PATH_SUFFIXES lib lib32 lib64)
if(ALBERTA_UTIL_LIB)
  set(_CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALBERTA_UTIL_LIB})
  check_library_exists(${ALBERTA_UTIL_LIB} alberta_calloc "" _ALBERTA_UTIL_LIB_FUNCTIONAL)
  if(NOT _ALBERTA_UTIL_LIB_FUNCTIONAL)
    message(WARN "Could not find symbol alberta_calloc in ${ALBERTA_UTIL_LIB}")
    cmake_pop_check_state()
    _dune_set_alberta(FALSE)
  endif(NOT _ALBERTA_UTIL_LIB_FUNCTIONAL)
else(ALBERTA_UTIL_LIB)
  message(WARNING "Could not find library alberta_util")
  cmake_pop_check_state()
  _dune_set_alberta(FALSE)
endif(ALBERTA_UTIL_LIB)

if(ALBERTA_LIBCHECK)
  foreach(dim RANGE 1 9)
    find_library(ALBERTA_${dim}D_LIB alberta_${dim}d PATH ${ALBERTA_DIR} PATH_SUFFIXES lib lib32 lib64
      Cache FILEPATH DOC "Alberta lib for ${dim}D" NO_DEFAULT_PATH)
    find_library(ALBERTA_${dim}D_LIB alberta_${dim}d  PATH_SUFFIXES lib lib32 lib64)
    if(ALBERTA_${dim}D_LIB)
      set(CMAKE_REQUIRED_LIBRARIES ${_CMAKE_REQUIRED_LIBRARIES_OLD} ${ALBERTA_${dim}D_LIB} ${ALBERTA_UTIL_LIB} ${DUNE_LIBS})
      check_library_exists(${ALBERTA_${dim}D_LIB} mesh_traverse "" ALBERTA_${dim}D_LIB_FOUND)
      if(ALBERTA_${dim}D_LIB_FOUND)
	list(APPEND ALBERTA_WORLD_DIMS ${dim})
      endif(ALBERTA_${dim}D_LIB_FOUND)
    endif(ALBERTA_${dim}D_LIB)
  endforeach(dim RANGE 1 9)
    message(STATUS "Found alberta libraries for dimensions ${ALBERTA_WORLD_DIMS}")
else(ALBERTA_LIBCHECK)
  message(WARNING "Disabled checking whether libalberta_Nd can be linked.")
endif(ALBERTA_LIBCHECK)
cmake_pop_check_state()

list(LENGTH ALBERTA_WORLD_DIMS _length)
if(length GREATER 0)
  _dune_set_alberta(TRUE)
endif(length GREATER 0)

if(ALBERTA_VERSION STREQUAL "2.0")
  set(DUNE_ALBERTA_VERSION 0x200)
elseif(ALBERTA_VERSION STREQUAL "3.0")
  set(DUNE_ALBERTA_VERSION 0x200)
else()
  message(FATAL_ERROR "Internal Inconsistency: Invalid Alberta version reported: $ALBERTA_VERSION.")
endif(ALBERTA_VERSION STREQUAL "2.0")

if(ALBERTA_FOUND)
  include(GridType)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM ASSERTION WORLDDIM == ALBERTA_DIM
    GRIDTYPE ALBERTAGRID DUNETYPE "Dune::AlbertaGrid< dimgrid >"
    HEADERS dune/grid/albertagrid.hh dune/grid/albertagrid/dgfparser.hh)
endif(ALBERTA_FOUND)

#set(AUTHOR_WARNING "Alberta test not yet finished. Disabling Alberta until test is complete")
#set(ALBERTA_FOUND ALBERTA_FOUND-NOTFOUND)
#set(HAVE_ALBERTA FALSE)
_dune_set_alberta(TRUE)

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
      message(WARNING "Alberta dimension not set. Please set it, e.g. use add_dune_alberta_flags(GRIDDIM 2 <target>)")
    else(NOT ADD_ALBERTA_WORLDDIM)
      list(FIND ALBERTA_WORLD_DIMS ${ADD_ALBERTA_WORLDDIM} _index)
      if(_index EQUAL -1)
	message(FATAL_ERROR "There is no alberta library for dimension ${ADD_ALBERTA_GRIDDIM}.")
      endif(_index EQUAL -1)
    endif(NOT ADD_ALBERTA_WORLDDIM)
    if(ADD_ALBERTA_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      set_property(DIRECTORY APPEND PROPERTY INCLUDE_DIRECTORIES ${ALBERTA_INCLUDES})
    else()
      if(NOT ADD_ALBERTA_OBJECT)
	# link to ALUGRID libraries
	foreach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
	  target_link_libraries(${_target} dunealbertagrid_${ADD_ALBERTA_GRIDDIM}d
	    ${ALBERTA_${ADD_ALBERTA_GRIDDIM}D_LIB}
	    dunegrid ${DUNE_LIBS} ${ALBERTA_UTIL_LIB} ${ALBERTA_EXTRA_LIBS})
	endforeach(_target ${ADD_ALBERTA_UNPARSED_ARGUMENTS})
      endif(NOT ADD_ALBERTA_OBJECT)
      set(_prefix TARGET)
      set_property(${_prefix} ${ADD_ALBERTA_UNPARSED_ARGUMENTS} APPEND
	PROPERTY INCLUDE_DIRECTORIES ${ALBERTA_INCLUDES})
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
