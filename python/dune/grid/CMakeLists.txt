add_python_targets(grid
  __init__
  _grid
  _grids
  grid_generator
  core
  map
  datahandle
)
dune_symlink_to_source_files(FILES _grid.cc)
