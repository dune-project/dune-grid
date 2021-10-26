import dune.grid

mshfile = "../../../doc/grids/gmsh/circle1storder.msh"
unstructuredGrid = dune.grid.ugGrid( (dune.grid.reader.gmsh, mshfile), dimgrid=2 )
