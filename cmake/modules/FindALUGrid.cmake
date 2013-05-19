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
#

set(ALUGRID_VERSION_REQUIRED 1.50)

# try to find ALUGrid's pkg-config file
if(NOT ALUGRID_ROOT)
  pkg_check_modules(PKG_ALUGRID "alugrid")
endif(NOT ALUGRID_ROOT)

# search and execute alugridversion, first only at positions given by the user
find_file(ALUGRID_VERSION_
  NAMES alugridversion
  PATHS
    ${ALUGRID_ROOT}
    ${PKG_ALUGRID_LIBRARY_DIRS}/..
    ${PKG_ALUGRID_INCLUDE_DIRS}/..
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH)
# search and execute alugridversion, including default paths
find_file(ALUGRID_VERSION
  NAMES alugridversion
  PATH_SUFFIXES bin)

# exit if we did not find alugridversion
if(NOT ALUGRID_VERSION)
  message("Could not find ALUGrid.")
  return()
endif(NOT ALUGRID_VERSION)

# we found alugrid, check if it is recent
# here we know we have ALUGRID_VERSION, otherwise we exit above
execute_process(COMMAND ${ALUGRID_VERSION} -c ${ALUGRID_VERSION_REQUIRED}
  OUTPUT_VARIABLE ALUGRID_VERSION)
if(ALUGRID_VERSION LESS 0)
  message(STATUS "ALUGrid version is less then ${ALUGRID_VERSION_REQUIRED}")
  _dune_set_alugrid(FALSE)
  return()
else(ALUGRID_VERSION LESS 0)
  message(STATUS "ALUGrid version is compatible")
endif(ALUGRID_VERSION LESS 0)

# looking for include path
find_path(ALUGRID_INCLUDE_PATH alugrid_serial.h
  PATHS
    ${ALUGRID_ROOT}
    ${PKG_ALUGRID_INCLUDE_DIRS}
  PATH_SUFFIXES
    include
    include/serial
  DOC "Include path of serial alugrid headers.")

find_library(ALUGRID_LIB alugrid
  PATHS
    ${ALUGRID_ROOT}
    ${PKG_ALUGRID_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib32 lib64
  DOC "ALUGrid library")

if(NOT ALUGRID_LIB)
  message("ALUGrid library not found.")
  return()
endif(NOT ALUGRID_LIB)

if(ALUGRID_INCLUDE_PATH)
  cmake_push_check_state() # store required flags
  set(ALUGRID_INCLUDES ${ALUGRID_INCLUDE_PATH} ${ALUGRID_INCLUDE_PATH}/serial
    ${ALUGRID_INCLUDE_PATH}/duneinterface)
  set(ALUGRID_DEFINITIONS -DENABLE_ALUGRID)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALUGRID_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALUGRID_LIB})
  set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} ${ALUGRID_DEFINITIONS})
  check_include_file_cxx(stlheaders.h HAVE_ALUGRID_SERIAL_H)
  if(HAVE_ALUGRID_SERIAL_H)
    check_cxx_source_compiles("
      #include <alugrid_defineparallel.h>
      #if ALU3DGRID_BUILD_FOR_PARALLEL == 0
      #error
      #endif
      int main(){}"
      ALUGRID_PARALLEL_FOUND)
  else(HAVE_ALUGRID_SERIAL_H)
    message("ALUGrid header not found")
  endif(HAVE_ALUGRID_SERIAL_H)
  cmake_pop_check_state()
else(ALUGRID_INCLUDE_PATH)
  message("alugrid_serial.h not found in ${ALUGRID_INCLUDE_PATH}.")
endif(ALUGRID_INCLUDE_PATH)

if(ALUGRID_PARALLEL_FOUND AND MPI_FOUND)
  # check for parallel ALUGrid
  find_path(ALUGRID_PARALLEL_INCLUDE_PATH alugrid_parallel.h
    PATHS ${ALUGRID_ROOT}
    PATH_SUFFIXES include include/parallel
    NO_DEFAULT_PATH)
  find_path(ALUGRID_PARALLEL_INCLUDE_PATH alugrid_parallel.h
    PATH_SUFFIXES include include/parallel)

  if(ALUGRID_PARALLEL_INCLUDE_PATH)
    cmake_push_check_state() # store rquired flags
    list(APPEND ALUGRID_INCLUDES ${ALUGRID_INCLUDE_PATH}/parallel)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ALUGRID_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ALUGRID_LIB})
    #set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} -DSOME_MORE_DEF)
    #check_include_file_cxx(alugrid_parallel.h ALUGRID_PARALLEL_HEADER)
    message(AUTHOR_WARNING "Include file alugrid_parallel.h is not checked because
the check produces linker errors.")

    set(HAVE_ALUGRID_PARALLEL_H 1
      CACHE INTERNAL "Have include alugrid_parallel.h ")
    cmake_pop_check_state()
    if(NOT HAVE_ALUGRID_PARALLEL_H)
      message("alugrid_parallel.h not found in ${ALUGRID_PARALLEL_INCLUDE_PATH}")
    endif(NOT HAVE_ALUGRID_PARALLEL_H)
  else(ALUGRID_PARALLEL_INCLUDE_PATH)
    message("alugrid_parallel.h not found.")
  endif(ALUGRID_PARALLEL_INCLUDE_PATH)
endif(ALUGRID_PARALLEL_FOUND AND MPI_FOUND)

check_library_exists(${ALUGRID_LIB} malloc "" ALULIB_FUNCTIONAL)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ALUGrid"
  DEFAULT_MSG
  ALULIB_FUNCTIONALY
  HAVE_ALUGRID_SERIAL_H)

mark_as_advanced(ALULIB_FUNCTIONALY)

# if both headers and library are found, store results

# finally set all variables
if(ALUGRID_FOUND)
  set(HAVE_ALUGRID TRUE)

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
endif(ALUGRID_FOUND)
