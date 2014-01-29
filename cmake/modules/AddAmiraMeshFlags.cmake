#
# Module providing convenience methods for compile binaries with AmiraMesh support.
#
# Provides the following functions:
#
# add_dune_amiramesh_flags(target1 target2 ...)
#
# adds AmiraMesh flags to the targets for compilation and linking
#

function(add_dune_amiramesh_flags _targets)
  if(AMIRAMESH_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} ${AMIRAMESH_LIBRARIES})
      get_target_property(_props ${_target} COMPILE_FLAGS)
      string(REPLACE "_props-NOTFOUND" "" _props "${_props}")
      set_target_properties(${_target} PROPERTIES COMPILE_FLAGS
        "${_props} ${AMIRAMESH_COMPILE_FLAGS}")
    endforeach(_target ${_targets})
  endif(AMIRAMESH_FOUND)
endfunction(add_dune_amiramesh_flags)
