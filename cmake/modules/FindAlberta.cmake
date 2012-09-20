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

macro(_dune_set_alberta val)
  set(ALBERTA_FOUND ${val})
  set(HAVE_ALBERTA $val})
endmacro(_dune_set_alberta val)

set(ALBERTA_LIBCHECK ON CACHE BOOL "Whether to try to link against libalberta_Nd.")
set(ALBERTA_DIR "" CACHE FILEPATH "Root directory of Alberta installation.")
set(ALBERTA_EXTRA_LIBS "ltdl;m" CACHE FILEPATH "Extra libraries needed by alberta for linking.")

# look for header alberta/alberta.h
find_path(ALBERTA_INCLUDE_DIR name alberta/alberta.h PATHS ${ALBERTA_DIR}
  PATH_SUFFIXES alberta include NO_DEFAULT_PATH
  DOC "Include path of Alberta")
find_path(ALBERTA_INCLUDE_DIR name alberta/alberta.h PATHS /usr/local /opt
  PATH_SUFFIXES alberta)

if(NOT ALBERTA_INCLUDE_DIR)
  _dune_set_alberta(FALSE)
  return()
endif(NOT ALBERTA_INCLUDE_DIR)

set(ALBERTA_INCLUDES ${ALBERTA_INCLUDE_DIR} ${ALBERTA_INCLUDE_DIR}/alberta)
set(ALBERTA_VERSION 2.0)

cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} -DDIM_OF_WORLD=3 -DDEL_INDEX=0)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALBERTA_INCLUDES})

check_include_files(alberta/alberta.h ALBERTA_FOUND)

if(NOT ALBERTA_FOUND)
  message(WARNING "Header alberta/alberta.h not found.")
  cmake_pop_check_state()
  _dune_set_alberta(FALSE)
  return()
endif(NOT ALBERTA_FOUND)

# check version, Alberta 2 or 3
include(CheckStructHasMember)
check_struct_has_member ("struct el_info" wall_bound alberta/alberta.h ALBERTA_IS_VERSION_3)
if(ALBERTA_IS_VERSION_3)
  set(ALBERTA_VERSION 3.0)
endif(ALBERTA_IS_VERSION_3)

# look for libraries
find_library(ALBERTA_UTIL_LIB
  NAMES alberta_util alberta_utilities
  PATHS ${ALBERTA_DIR}
  PATH_SUFFIXES lib lib32 lib64
  NO_DEFAULT_PATH)
find_library(ALBERTA_UTIL_LIB
  NAMES alberta_util alberta_utilities
  PATH_SUFFIXES lib lib32 lib64)

if(ALBERTA_UTIL_LIB)
  set(_CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALBERTA_UTIL_LIB})
  check_library_exists(${ALBERTA_UTIL_LIB} alberta_calloc "" _ALBERTA_UTIL_LIB_FUNCTIONAL)
  if(NOT _ALBERTA_UTIL_LIB_FUNCTIONAL)
    message(WARNING "Could not find symbol alberta_calloc in ${ALBERTA_UTIL_LIB}")
    cmake_pop_check_state()
    _dune_set_alberta(FALSE)
    return()
  endif(NOT _ALBERTA_UTIL_LIB_FUNCTIONAL)
else(ALBERTA_UTIL_LIB)
  message(WARNING "Could not find library alberta_util or alberta_utilities")
  cmake_pop_check_state()
  _dune_set_alberta(FALSE)
endif(ALBERTA_UTIL_LIB)

set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALBERTA_UTIL_LIB} ${ALBERTA_EXTRA_LIBS})

# check for which dimensions are supported Alberta installation
if(ALBERTA_LIBCHECK)
  foreach(dim RANGE 1 9)
    find_library(ALBERTA_${dim}D_LIB alberta_${dim}d
      PATHS ${ALBERTA_DIR}
      PATH_SUFFIXES lib lib32 lib64
      Cache FILEPATH DOC "Alberta lib for ${dim}D" NO_DEFAULT_PATH)
    find_library(ALBERTA_${dim}D_LIB alberta_${dim}d  PATH_SUFFIXES lib lib32 lib64)
    if(ALBERTA_${dim}D_LIB)
      #set(CMAKE_REQUIRED_LIBRARIES ${_CMAKE_REQUIRED_LIBRARIES_OLD} ${ALBERTA_${dim}D_LIB} ${ALBERTA_UTIL_LIB} ${DUNE_LIBS})
      check_library_exists(${ALBERTA_${dim}D_LIB} mesh_traverse "" ALBERTA_${dim}D_LIB_FOUND)
      if(ALBERTA_${dim}D_LIB_FOUND)
        list(APPEND ALBERTA_WORLD_DIMS ${dim})
      endif(ALBERTA_${dim}D_LIB_FOUND)
    endif(ALBERTA_${dim}D_LIB)
  endforeach(dim RANGE 1 9)
  message(STATUS "Found alberta libraries for dimensions ${ALBERTA_WORLD_DIMS}")
endif(ALBERTA_LIBCHECK)
cmake_pop_check_state()

list(LENGTH ALBERTA_WORLD_DIMS _length)
if(length GREATER 0)
  _dune_set_alberta(TRUE)
endif(length GREATER 0)

if(ALBERTA_VERSION STREQUAL "2.0")
  set(DUNE_ALBERTA_VERSION 0x200)
elseif(ALBERTA_VERSION STREQUAL "3.0")
  set(DUNE_ALBERTA_VERSION 0x300)
else()
  message(WARNING "Internal Inconsistency: Invalid Alberta version reported: $ALBERTA_VERSION.")
  cmake_pop_check_state()
  _dune_set_alberta(FALSE)
endif(ALBERTA_VERSION STREQUAL "2.0")

if(ALBERTA_FOUND)
  include(GridType)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM ASSERTION WORLDDIM == ALBERTA_DIM
    GRIDTYPE ALBERTAGRID DUNETYPE "Dune::AlbertaGrid< dimgrid >"
    HEADERS dune/grid/albertagrid.hh dune/grid/albertagrid/dgfparser.hh)
  _dune_set_alberta(TRUE)
endif(ALBERTA_FOUND)
