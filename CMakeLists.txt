add_python_targets(grid
  __init__
  _grids
  grid_generator
  core
  map
  datahandle
)

dune_add_pybind11_module(NAME _grid)
