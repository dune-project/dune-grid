
macro(_dune_set_alberta val)
  set(ALBERTA_FOUND ${val})
  set(HAVE_ALBERTA ${val})
endmacro(_dune_set_alberta val)

set(ALBERTA_LIBCHECK ON CACHE BOOL "Whether to try to link against libalberta_Nd.")
set(ALBERTA_ROOT "" CACHE FILEPATH "Root directory of Alberta installation.")
set(ALBERTA_EXTRA_LIBS "ltdl;m" CACHE FILEPATH "Extra libraries needed by alberta for linking.")

# look for header alberta/alberta.h
find_path(ALBERTA_INCLUDE_DIR
  NAMES alberta/alberta.h
  PATHS ${ALBERTA_ROOT}
  PATH_SUFFIXES alberta include NO_DEFAULT_PATH
  DOC "Include path of Alberta")
find_path(ALBERTA_INCLUDE_DIR
  NAMES
  alberta/alberta.h
  PATHS /usr/local /opt
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
  PATHS ${ALBERTA_ROOT}
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
      PATHS ${ALBERTA_ROOT}
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

#add all alberta grid related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_alberta_flags
if(ALBERTA_FOUND)
  foreach(dir ${ALBERTA_INCLUDES})
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()
