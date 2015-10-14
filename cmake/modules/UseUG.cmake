# This module first tests for UG and then sets the necessary flags
# and config.h defines. If UG is found UG_FOUND will be true.
#
# .. cmake_module::
#
#    Checks for presence and usability of the UG library.
#
#    You may set the following variables to configure this module:
#
#    :ref:`UG_ROOT`
#       Path list to search for UG
#
#    This module sets the following variables:
#
#    :code:`UG_FOUND`
#       True if the UG library was found
#
#    :code:`UG_VERSION`
#       The version string of the found UG library
#
#    :code:`UG_INCLUDES`
#       The include directories needed or UG
#
#    :code:`UG_LIBRARIES`
#       The libraries, UG needs to link to
#
#    .. note::
#       This module is not called `FindUG`, because UG does ship the module
#       of that name and we would otherwise shadow that one.
#
# .. cmake_variable:: UG_ROOT
#
#    You may set this variable to have :ref:`UseUG` look
#    for the UG package in the given path before inspecting
#    system paths.
#
# .. cmake_function:: add_dune_ug_flags
#
#    .. cmake_param:: targets
#       :single:
#       :positional:
#       :required:
#
#       The targets to add the UG flags to.
#
#    .. cmake_param:: SOURCE_ONLY
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are source files.
#
#    .. cmake_param:: OBJECT
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are object libraries.
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

find_package(UG 3.11.0
  NO_MODULE QUIET
  NO_DEFAULT_PATH)
find_package(UG 3.11.0
  NO_MODULE
  NO_SYSTEM_ENVIRONMENT_PATH)

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

  # register all UG related flags
  set(UG_DEFINITIONS "ENABLE_UG=1")
  if(UG_PARALLEL STREQUAL "yes")
    set(UG_DEFINITIONS "ENABLE_UG=1;ModelP")
  endif()
  dune_register_package_flags(COMPILE_DEFINITIONS "${UG_DEFINITIONS}"
                              INCLUDE_DIRS "${UG_INCLUDES}"
                              LIBRARIES "dunegrid;${UG_LIBRARIES};${DUNE_LIBS}")

  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of UG ${UG_VERSION} succeeded:\n"
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
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT;NO_LINK_DUNEGRID" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
    else()
      if(NOT ADD_UG_OBJECT)
        foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
          if(NOT ADD_UG_NO_LINK_DUNEGRID)
            target_link_libraries(${_target} dunegrid)
          endif()
          target_link_libraries(${_target}
            ${UG_LIBRARIES} ${DUNE_LIBS})
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
      if(NOT ADD_UG_NO_LINK_DUNEGRID)
        set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
          APPEND PROPERTY
          LINK_LIBRARIES dunegrid)
      endif()
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        LINK_LIBRARIES ${UG_LIBRARIES} ${DUNE_LIBS})
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
