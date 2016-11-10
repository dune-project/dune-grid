# This module first tests for the old UG and then sets the necessary flags
# and config.h defines. If UG is found UG_FOUND will be true.
# If dune-uggrid was found it only adds the dgf magic to config.h
# and makes add_dune_ug_flags available to stay compatible.
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

# Support for external UG is deprecated in Dune 2.5. We will keep it in Dune 2.6-git as long as feasible.
set_package_info("UG" "External UG grid, superseded by dune-uggrid" "http://www.iwr.uni-heidelberg.de/frame/iwrwikiequipment/software/ug")

if(NOT dune-uggrid_FOUND)
  if(UG_ROOT AND NOT UG_DIR)
    # define the directory where the config file resides
    if(EXISTS "${UG_ROOT}/lib/cmake/ug/ug-config.cmake")
      set(UG_DIR ${UG_ROOT}/lib/cmake/ug)
    elseif(EXISTS "${UG_ROOT}/lib64/cmake/ug/ug-config.cmake")
      set(UG_DIR ${UG_ROOT}/lib64/cmake/ug)
    else()
      message(WARNING "Could not find file ug-config.cmake relative to given UG_ROOT")
    endif()
  endif()

  find_package(UG 3.11.0
    NO_MODULE QUIET
    NO_DEFAULT_PATH)
  find_package(UG 3.11.0
    NO_MODULE
    NO_SYSTEM_ENVIRONMENT_PATH)

  set(HAVE_UG ${UG_FOUND})

  if(UG_FOUND)
    # parse patch level: last number in UG version string is DUNE patch level
    string(REGEX MATCH "[0-9]*$" UG_DUNE_PATCHLEVEL ${UG_VERSION})

    if (UG_VERSION VERSION_GREATER 3.13.0
        OR UG_VERSION VERSION_EQUAL 3.13.0)
      list(APPEND UG_DEFINITIONS "UG_USE_NEW_DIMENSION_DEFINES")
    endif()

    dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID ASSERTION GRIDDIM == WORLDDIM
        DUNETYPE "Dune::UGGrid< dimgrid >"
        HEADERS dune/grid/uggrid.hh dune/grid/io/file/dgfparser/dgfug.hh)

    # Overwrite flags by hand
    set(UG_LIBRARIES "")
    set(_paths "${prefix}")

    # Find out the full path to the libs
    foreach(entry ${UG_LIBRARY_FLAGS} -L/bla)
      string(REGEX REPLACE "^-L([a-zA-Z/-_]+)" "\\1" _path ${entry})
      list(APPEND _paths ${_path})
    endforeach()

    foreach(lib ugS2 ugS3 devS)
        set(full_path "full_path-NOTFOUND")
        find_library(full_path ${lib} PATHS ${_paths} NO_DEFAULT_PATH)
        if(full_path)
          list(APPEND UG_LIBRARIES ${full_path})
        endif()
    endforeach()

    # register all UG related flags
    list(APPEND UG_DEFINITIONS "ENABLE_UG=1")
    if(UG_PARALLEL STREQUAL "yes")
      list(APPEND UG_DEFINITIONS "ModelP")
    endif()

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
endif()

# Add dgf magic to config.h and register flags
if(UG_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "${UG_DEFINITIONS}"
                              INCLUDE_DIRS "${UG_INCLUDES}"
                              LIBRARIES "dunegrid;${UG_LIBRARIES};${DUNE_LIBS}")
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID ASSERTION GRIDDIM == WORLDDIM
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS dune/grid/uggrid.hh dune/grid/io/file/dgfparser/dgfug.hh)

    # Support for external UG is deprecated in Dune 2.5. We will keep it in Dune 2.6-git as long as feasible.
    message(WARNING "The support of UG as an external library is deprecated in Dune 2.5. Use dune-uggrid instead.")
endif()

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
        endforeach()
      endif()
      set(_prefix TARGET)
    endif()

    # Add compiler flags
    include_directories(${UG_INCLUDES})
    set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
      APPEND PROPERTY
      COMPILE_DEFINITIONS ENABLE_UG)
    if (UG_VERSION VERSION_GREATER 3.13.0
        OR UG_VERSION VERSION_EQUAL 3.13.0)
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS}
        APPEND PROPERTY
        COMPILE_DEFINITIONS UG_USE_NEW_DIMENSION_DEFINES)
    endif()
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
    endif()
  endif()
endfunction(add_dune_ug_flags)
