## -*- autoconf -*-
AC_DEFUN([DUNE_GRID_CHECKS],[
  AC_REQUIRE([DUNE_GRID_DIMENSION])
  AC_REQUIRE([DUNE_PATH_GRAPE])
  AC_REQUIRE([DUNE_PATH_PARMETIS])
  AC_REQUIRE([DUNE_PATH_ALBERTA])
  AC_REQUIRE([DUNE_PATH_UG])
  AC_REQUIRE([DUNE_PATH_AMIRAMESH])
  AC_REQUIRE([DUNE_PATH_PSURFACE])
  AC_REQUIRE([DUNE_PATH_ALUGRID])
  AC_REQUIRE([DUNE_EXPERIMENTAL_GRID_EXTENSIONS])

  DUNE_DEFINE_GRIDTYPE([ONEDGRID],[(GRIDDIM == 1) && (WORLDDIM == 1)],[Dune::OneDGrid],[dune/grid/onedgrid.hh],[dune/grid/io/file/dgfparser/dgfoned.hh])
  DUNE_DEFINE_GRIDTYPE([SGRID],[],[Dune::SGrid< dimgrid, dimworld >],[dune/grid/sgrid.hh],[dune/grid/io/file/dgfparser/dgfs.hh])
  DUNE_DEFINE_GRIDTYPE([YASPGRID],[GRIDDIM == WORLDDIM],[Dune::YaspGrid< dimgrid >],[dune/grid/yaspgrid.hh],[dune/grid/io/file/dgfparser/dgfyasp.hh])
])

AC_DEFUN([DUNE_GRID_CHECK_MODULE],[
  DUNE_CHECK_MODULES([dune-grid], [grid/onedgrid.hh],[dnl
  std::vector<Dune::OneDGrid::ctype> coords;
  Dune::OneDGrid grid(coords);
  return grid.lbegin<0>(0) == grid.lend<0>(0);])
])
