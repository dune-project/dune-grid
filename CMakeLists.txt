add_python_targets(grid
  __init__
  _create
  _grids
  grid_generator
  core
  map
)

add_library(_grid SHARED _grid.cc)
set_target_properties(_grid PROPERTIES PREFIX "")
