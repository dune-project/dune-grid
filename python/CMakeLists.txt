add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  DEPENDS _grid
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE HAVE_ALBERTA HAVE_UG
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
