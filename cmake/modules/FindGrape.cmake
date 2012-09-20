# Module that checks whether Grape is available
#
# Accepts the following input variable
# GRAPE_PREFIX: Prefix under which Grape is installed
# GRAPE_INCLUDE_DIRS: Include directories for Grape
# GRAPE_LIBRARY: Full path to Grape library
#
# The following variable will be set:
# GRAPE_FOUND: whether Grape is available
# GRAPE_INCLUDE_DIRS: Include directories for Grape
# GRAPE_LIBRARIES: Full path to libraries needed to link
#   to Grape
#
# Provides the function
# add_dune_grape_flags( [OBJECT | SOURCE_ONLY] target1 ...)
#   that sets all necessary flags needed for compilation and linking.
#


message(AUTHOR_WARNING "TODO: Test Grape test and macros")
set(GRAPE_FOUND GRAPE_FOUND-NOTFOUND)

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

find_package(X11)
find_package(OpenGL)
set(GRAPE_PREFIX "/usr/local/grape/" CACHE FILEPATH "Prefix directory, where GRAPE is installed")

if(X11_FOUND)
   # find header in user supplied directory
  find_path(GRAPE_INCLUDE_DIRS grape.h PATHS ${GRAPE_PREFIX}
    DOC "Include directory with Grape header files" NO_DEFAULT_PATH)
  find_path(GRAPE_INCLUDE_DIRS grape.h) #standard directories

  # check header usability
  include(CMakePushCheckState)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -DENABLE_GRAPE")
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${OPENGL_INCLUDE_DIR} ${GRAPE_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${OPENGL_LIBRARIES} ${XEXT_LIB} ${CMAKE_REQUIRED_LIBRARIES})
  check_include_files(grape.h _GRAPE_HEADER_USABLE)

  # find library
  find_library(GRAPE_LIBRARY NAMES grape
    PATHS ${GRAPE_PREFIX} NO_DEFAULT_PATH
    DOC "Full path to grape library.")
  find_library(GRAPE_LIBRARY NAMES grape)

  find_path(GRAPE_LIBRARY_PATH GRAPE_LIBRARY NO_DEFAULT_PATH)

  include(CheckLibraryExists)
  check_library_exists(gr grape GRAPE_PREFIX _GRAPE_LIB_FUNCTIONAL)
  cmake_pop_check_state()

  if(_GRAPE_LIB_FUNCTIONAL)
    set(GRAPE_LIBRARIES ${GRAPE_LIBRARY} ${OPENGL_LIBRARIES} ${XEXT_LIB})
  endif(_GRAPE_LIB_FUNCTIONAL)
else(X11_FOUND)
  message(WARNING "[X libraries were not found and therefore not Grape check possible!")
endif(X11_FOUND)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "GRAPE"
  DEFAULT_MSG
  GRAPE_INCLUDE_DIRS
  GRAPE_LIBRARIES
  _GRAPE_LIB_FUNCTIONAL
  _GRAPE_HEADER_USABLE
)
set(HAVE_GRAPE ${GRAPE_FOUND})
mark_as_advanced(GRAPE_INCLUDE_DIRS GRAPE_LIBRARIES)
