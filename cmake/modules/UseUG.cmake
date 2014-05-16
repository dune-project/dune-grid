#
# This module first tests for UG and then sets the necessary flags
# and config.h defines. If UG is found UG_FOUND will be true.
#
# Provides the following function:
#
# add_dune_ug_flags(<target1> [<target2> ...] [SOURCE_ONLY] [OBJECT])
#
# This functions sets the necessary flags for a program that uses UG.
# The option SOURCE_ONLY indicates that the targets are source files.
# The option OBJECT indicates that the targets are object libraries.
#
if(UG_ROOT AND NOT UG_DIR)
  # define the directory where the config file resides
  if(EXISTS "${UG_ROOT}/lib/cmake/ug/ug-config.cmake")
    set(UG_DIR ${UG_ROOT}/lib/cmake/ug)
  elseif(EXISTS "${UG_ROOT}/lib64/cmake/ug/ug-config.cmake")
    set(UG_DIR ${UG_ROOT}/lib64/cmake/ug)
  else()
    message(WARNING "Could not find file ug-config.cmake relative to given UG_ROOT")
  endif()
endif(UG_ROOT AND NOT UG_DIR)

find_package(UG 3.9.1
  NO_MODULE QUIET
  NO_DEFAULT_PATH)
find_package(UG 3.9.1
  NO_MODULE)

if(UG_FOUND AND (NOT UG_FOR_DUNE STREQUAL "yes"))
  message(WARNING "UG was not configured for DUNE. Did you pass --enable-dune to its configure?")
  set(UG_FOUND "UG_FOUND-NOTFOUND")
endif()

set(HAVE_UG ${UG_FOUND})

if(UG_FOUND)
  # parse patch level: last number in UG version string is DUNE patch level
  string(REGEX MATCH "[0-9]*$" UG_DUNE_PATCHLEVEL ${UG_VERSION})

  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID ASSERTION GRIDDIM == WORLDDIM
      DUNETYPE "Dune::UGGrid< dimgrid >"
      HEADERS dune/grid/uggrid.hh dune/grid/io/file/dgfparser/dgfug.hh)

  # Remove the following as soon as we absolutely require patch10 or higher
  if(${UG_DUNE_PATCHLEVEL} GREATER 9)
    set(HAVE_UG_PATCH10 1)
  else()
    set(HAVE_UG_PATCH10 0)
  endif(${UG_DUNE_PATCHLEVEL} GREATER 9)

  #Overwrite flags by hand (like for autoconf).
  set(UG_LIBRARIES)
  set(paths "${prefix}")

  #Find out the full path to the libs.
  foreach(entry ${UG_LIBRARY_FLAGS} -L/bla)
    string(REGEX REPLACE "^-L([a-zA-Z/-_]+)" "\\1" _path ${entry})
    list(APPEND _paths ${_path})
  endforeach(entry {UG_LIBRARY_FLAGS})

  foreach(lib ugS2 ugS3 devS)
      set(full_path "full_path-NOTFOUND")
      find_library(full_path ${lib} PATHS ${_paths} NO_DEFAULT_PATH)
      if(full_path)
        list(APPEND UG_LIBRARIES ${full_path})
      endif(full_path)
  endforeach(lib ugS2 ugS3 devS)

  # add all UG related flags to ALL_PKG_FLAGS, this must happen
  # regardless of a target using add_dune_ug_flags
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-DENABLE_UG")
  foreach(dir ${UG_INCLUDES})
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of UG ${UG_VERSION} succeded:\n"
    "Include directories: ${UG_INCLUDES}\n"
    "Libraries: ${UG_LIBRARIES}\n\n")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of UG failed:\n"
    "Include directories: ${UG_INCLUDES}\n"
    "Libraries: ${UG_LIBRARIES}\n\n")
endif()

# output whether UG found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "UG"
  DEFAULT_MSG
  UG_DIR
  HAVE_UG
)

# Add flags to targets
function(add_dune_ug_flags)
  if(UG_FOUND)
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
    else()
      if(NOT ADD_UG_OBJECT)
        foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target}
            dunegrid ${UG_LIBRARIES})
        endforeach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
      endif()
      set(_prefix TARGET)
    endif()

    # Add compiler flags
    include_directories(${UG_INCLUDES})
    set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
      APPEND PROPERTY
      COMPILE_DEFINITIONS ENABLE_UG)
    # Add linker arguments
    if(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        LINK_LIBRARIES ${UG_LIBRARIES})
    endif(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
    if(UG_PARALLEL STREQUAL "yes")
      # Add modelp
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        COMPILE_DEFINITIONS ModelP)
      # Add mpi flags.
      add_dune_mpi_flags(${ADD_UG_UNPARSED_ARGUMENTS} ${_source_only})
    endif(UG_PARALLEL STREQUAL "yes")
  endif(UG_FOUND)
endfunction(add_dune_ug_flags)
