# Module providing convenience methods for compile binaries with psurface support.
#
# .. cmake_function:: add_dune_psurface_flags
#
#    .. cmake_param:: targets
#       :single:
#       :required:
#       :positional:
#
#       the targets to add the Psurface flags to.
#

function(add_dune_psurface_flags _targets)
  if(PSURFACE_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} ${PSURFACE_DUNE_LIBRARIES})
      get_target_property(_props ${_target} COMPILE_FLAGS)
      string(REPLACE "_props-NOTFOUND" "" _props "${_props}")
      set_target_properties(${_target} PROPERTIES COMPILE_FLAGS
        "${_props} ${PSURFACE_DUNE_COMPILE_FLAGS}")
    endforeach(_target ${_targets})
  endif(PSURFACE_FOUND)
endfunction(add_dune_psurface_flags)
