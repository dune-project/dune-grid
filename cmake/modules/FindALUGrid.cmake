#
# Module that checks whether ALUGrid is available.
# ALUGrid must be at least version 1.50.
#
# Variables used by this module which you may want to set:
# ALUGRID_ROOT   Path list to search for ALUGrid
#
# Sets the follwing variable:
#
# ALUGRID_FOUND           True if ALUGrid available.
# HAVE_ALUGRID            True if ALUGrid available.
# HAVE_ALUGRID_SERIAL_H   1 if serial header found.
# HAVE_ALUGRID_PARALLEL_H 1 if parallel header found, too.
#

set(ALUGRID_VERSION_REQUIRED 1.50)

# try to find ALUGrid's pkg-config file
if(NOT ALUGRID_ROOT)
  pkg_check_modules(PKG_ALUGRID "alugrid")
endif(NOT ALUGRID_ROOT)

# try manually, if ALUGrid's pkg-config file not found
if(NOT PKG_ALUGRID_FOUND)
  # first only at positions given by the user
  find_path(PATH_PKG_ALUGRID
    NAMES "alugrid.pc"
    PATHS ${ALUGRID_ROOT}
    PATH_SUFFIXES lib/pkgconfig lib32/pkgconfig lib64/pkgconfig
    NO_DEFAULT_PATH)
  # including default paths
  find_path(PATH_PKG_ALUGRID
    NAMES "alugrid.pc"
    PATH_SUFFIXES lib/pkgconfig lib32/pkgconfig lib64/pkgconfig)

  # try again with path temporarilly added to PKG_CONFIG_PATH
  set(REM_PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
  set(ENV{PKG_CONFIG_PATH} "${ALUGRID_ROOT}:${PATH_PKG_ALUGRID}:${PKG_CONFIG_PATH}")
  pkg_check_modules(PKG_ALUGRID "alugrid")
  set(ENV{PKG_CONFIG_PATH} REM_PKG_CONFIG_PATH)
endif(NOT PKG_ALUGRID_FOUND)

# check whether ALUGrid version is recent enough
if(PKG_ALUGRID_FOUND)
  if(ALUGRID_VERSION_REQUIRED VERSION_GREATER "${PKG_ALUGRID_VERSION}")
    message(STATUS "Found ALUGrid ${PKG_ALUGRID_VERSION}, but version ${ALUGRID_VERSION_REQUIRED} is required")
  else()
    set(ALUGRID_VERSION "${PKG_ALUGRID_VERSION}")
  endif()
endif(PKG_ALUGRID_FOUND)

# look for include path
if(PKG_ALUGRID_FOUND)
  find_path(ALUGRID_INCLUDE_PATH alugrid_serial.h
    PATHS
      ${ALUGRID_ROOT}
      ${PKG_ALUGRID_INCLUDE_DIRS}
    PATH_SUFFIXES
      include
      include/serial
    DOC "Include path of serial alugrid headers.")
endif(PKG_ALUGRID_FOUND)

# look for library path
if(ALUGRID_INCLUDE_PATH)
  find_library(ALUGRID_LIB alugrid
    PATHS
      ${ALUGRID_ROOT}
      ${PKG_ALUGRID_LIBRARY_DIRS}
    PATH_SUFFIXES lib lib32 lib64
    DOC "ALUGrid library")
endif(ALUGRID_INCLUDE_PATH)

# check whether library works
if(ALUGRID_LIB)
  get_filename_component(ALUGRID_LIB_PATH ${ALUGRID_LIB} PATH)
  check_library_exists(alugrid malloc ${ALUGRID_LIB_PATH} ALULIB_FUNCTIONAL)
endif(ALUGRID_LIB)

set(ALUGRID_INCLUDES ${ALUGRID_INCLUDE_PATH} ${ALUGRID_INCLUDE_PATH}/serial
  ${ALUGRID_INCLUDE_PATH}/duneinterface)
set(ALUGRID_DEFINITIONS -DENABLE_ALUGRID)

include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALUGRID_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALUGRID_LIB})
set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} ${ALUGRID_DEFINITIONS})

# try to use a header
if(ALUGRID_LIB)
  check_include_file_cxx(stlheaders.h HAVE_ALUGRID_SERIAL_H)
endif(ALUGRID_LIB)

# check whether it is parallel ALUGrid
if(HAVE_ALUGRID_SERIAL_H AND MPI_CXX_FOUND)
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles("
    #include <alugrid_defineparallel.h>
    #if ALU3DGRID_BUILD_FOR_PARALLEL == 0
    #error
    #endif
    int main(){}"
    ALUGRID_PARALLEL_FOUND)
endif(HAVE_ALUGRID_SERIAL_H AND MPI_CXX_FOUND)

# try to use parallel ALUGrid
if(ALUGRID_PARALLEL_FOUND AND MPI_CXX_FOUND)
  # find path to parallel headers
  find_path(ALUGRID_PARALLEL_INCLUDE_PATH alugrid_parallel.h
    PATHS
      ${ALUGRID_ROOT}
      ${PKG_ALUGRID_INCLUDE_DIRS}
    PATH_SUFFIXES include include/parallel
    NO_DEFAULT_PATH)
  find_path(ALUGRID_PARALLEL_INCLUDE_PATH alugrid_parallel.h
    PATH_SUFFIXES include include/parallel)

  if(ALUGRID_PARALLEL_INCLUDE_PATH)
    list(APPEND ALUGRID_INCLUDES ${ALUGRID_INCLUDE_PATH}/parallel)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALUGRID_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_DUNE_LIBRARIES})
    check_cxx_source_compiles("
      #include <alugrid_parallel.h>
      int main(){}"
      ALUGRID_PARALLEL_WORKS)

    if(ALUGRID_PARALLEL_WORKS)
      set(HAVE_ALUGRID_PARALLEL_H 1
        CACHE INTERNAL "Have include alugrid_parallel.h")
    endif(ALUGRID_PARALLEL_WORKS)
  endif(ALUGRID_PARALLEL_INCLUDE_PATH)
endif(ALUGRID_PARALLEL_FOUND AND MPI_CXX_FOUND)

cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ALUGrid"
  DEFAULT_MSG
  ALUGRID_VERSION
  ALUGRID_INCLUDE_PATH ALUGRID_LIB ALUGRID_LIB_PATH ALULIB_FUNCTIONAL
  HAVE_ALUGRID_SERIAL_H)

mark_as_advanced(PATH_PKG_ALUGRID REM_PKG_CONFIG_PATH
  ALUGRID_INCLUDE_PATH ALUGRID_LIB ALUGRID_LIB_PATH ALULIB_FUNCTIONAL
  ALUGRID_PARALLEL_INCLUDE_PATH ALUGRID_PARALLEL_WORKS)

# set HAVE_ALUGRID for config.h
set(HAVE_ALUGRID ${ALUGRID_FOUND})

# finally set all variables
if(ALUGRID_FOUND)
  set(ALUGRID_LIBRARIES    ${ALUGRID_LIB})

  include(GridType)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_CONFORM
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, simplex, conforming >"
    HEADERS dune/grid/alugrid.hh dune/grid/io/file/dgfparser/dgfalu.hh)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_CUBE
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, cube, nonconforming >"
    HEADERS dune/grid/alugrid.hh dune/grid/io/file/dgfparser/dgfalu.hh)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_SIMPLEX
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, simplex, nonconforming >"
    HEADERS dune/grid/alugrid.hh dune/grid/io/file/dgfparser/dgfalu.hh)

  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ALUGrid ${ALUGRID_VERSION} succeded:\n"
    "Include directory: ${ALUGRID_INCLUDES}\n"
    "Library directory: ${ALUGRID_LIBRARIES}\n\n")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ALUGrid failed.\n\n")
endif(ALUGRID_FOUND)

#add all alugrid related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_alugrid_flags
if(ALUGRID_FOUND)
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_DEFS "ENABLE_ALUGRID=1")
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_INCS "${ALUGRID_INCLUDES}")
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_LIBS "${ALUGRID_LIBRARIES}")
endif()
