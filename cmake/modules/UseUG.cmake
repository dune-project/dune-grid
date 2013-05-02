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

find_package(UG 3.9.1)
if(NOT CMAKE_DISABLE_FIND_PACKAGE_UG)
  if(NOT UG_FOUND)
    message(WARNING "CMake will only find the most current UG. Maybe you need to upgrade?")
  endif(NOT UG_FOUND)
  message(AUTHOR_WARNING "We need to test for the patch level, too!")
endif(NOT CMAKE_DISABLE_FIND_PACKAGE_UG)
set(HAVE_UG ${UG_FOUND})

dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID ASSERTION GRIDDIM == WORLDDIM
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS dune/grid/uggrid.hh dune/grid/io/file/dgfparser/dgfug.hh)

#Overwrite flags by hand (like for autoconf).
set(UG_LIBRARIES dunegrid)
set(paths "/home/mblatt/src/dune/3rdParty/grids/ug-install")

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

function(add_dune_ug_flags)
  if(UG_FOUND)
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${UG_INCLUDES})
    else()
      set(_prefix TARGET)
      if(ADD_UG_OBJECT)
	set(_prefix TARGET)
      else(ADD_UG_OBJECT)
	foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
	  target_link_libraries(${_target} ${UG_LIBRARIES} ${DUNE_LIBS})
	endforeach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
      endif(ADD_UG_OBJECT)
      include_directories(${UG_INCLUDES})
    endif()

    set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND PROPERTY COMPILE_DEFINITIONS ENABLE_UG)
    if(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND PROPERTY LINK_LIBRARIES ${UG_LIBRARIES} dunegrid ${DUNE_LIBS})
    endif(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
  endif(UG_FOUND)
endfunction(add_dune_ug_flags)
