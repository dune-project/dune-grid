add_python_targets(grid
  __init__
  _create
  grid_generator
  yasp
  core
  map
  alberta
  oned
  sp
  ug
)

add_library(_grid SHARED _grid.cc)
set_target_properties(_grid PROPERTIES PREFIX "")
