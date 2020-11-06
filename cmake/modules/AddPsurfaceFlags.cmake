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

# set HAVE_PSURFACE for config.h
set(HAVE_PSURFACE ${Psurface_FOUND})

# register all psurface related flags
if(Psurface_FOUND)
  dune_register_package_flags(LIBRARIES "Psurface::Psurface")
endif()

function(add_dune_psurface_flags _targets)
  if(Psurface_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Psurface::Psurface)
    endforeach(_target ${_targets})
  endif(Psurface_FOUND)
endfunction(add_dune_psurface_flags)
