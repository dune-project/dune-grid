// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*
 * Check leaf intersection geometry on UGGrid.
 *
 * Reference: https://gitlab.dune-project.org/core/dune-grid/issues/53
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk.hh>
#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/grid.hh>
#endif

const double EPSILON = 1e-12;

using Dune::TestSuite;

template<class Grid>
TestSuite testGrid(const Grid& grid, const std::string& name) {
  using std::sqrt;
  TestSuite test;

  // now iterate over intersections
  const auto& gv = leafGridView(grid);
  for (const auto& e : elements(gv)) {
    for (const auto& is : intersections(gv, e)) {
      if (is.boundary())
        continue;

      // compute the edge lengths of the inner and outer element
      auto insideLength = sqrt(is.inside().geometry().volume());
      auto outsideLength = sqrt(is.outside().geometry().volume());
      // length of intersection
      auto intersectionLength = is.geometry().volume();

      // intersection should not be longer than any full edge
      const bool check = intersectionLength < insideLength+EPSILON && intersectionLength < outsideLength+EPSILON;
      test.check(check)
        << "=== " << name << ": Intersection is longer than actual edge! ===\n"
        << "  inside center: " << is.inside().geometry().center() << "\n"
        << " outside center: " << is.outside().geometry().center() << "\n"
        << " lengths: intersection: " << intersectionLength << " inside: " << insideLength << " outside: " << outsideLength << "\n";
    }
  }

  // output mesh as VTK (uncomment if needed)
  // auto vtk = Dune::VTKWriter<typename Grid::LeafGridView>(gv);
  // vtk.write(s);

  return test;
}

template<class Grid>
TestSuite testRefinement1(Grid& grid, const std::string& name)
{
  const auto& gv = leafGridView(grid);

  // Make some pseudo-adaptive refinements
  for (int i = 0; i < 5; ++i)
  {
    for (auto e : elements(gv)) {
      if (gv.indexSet().subIndex(e,0,0) % 5 == 0) // mark some elements for refinement
        grid.mark(1, e);
    }
    grid.adapt();
    grid.postAdapt();
  }

  return testGrid(grid, name);
}

template<typename Grid>
TestSuite testRefinement2(Grid& grid, const std::string& name)
{
  using Coordinate = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;

  // List of list of elements to refine in each refinement cycle.
  // The elements are specified using their cell center.
  const std::vector< std::vector<Coordinate> > refine = {
    {{0.5, 0.5}},
    {{0.75, 0.25}},
    {{0.875, 0.375}},
    {{0.9375, 0.4375}, {0.875, 0.125}}
  };

  const auto& gv = leafGridView(grid);

  for (const auto& xs : refine) {
    for (const auto& e : elements(gv)) {
      for (const auto& x : xs) {
        auto d = e.geometry().center() - x;
        if (d.two_norm() < EPSILON)
          grid.mark(1, e);
      }
    }
    grid.adapt();
    grid.postAdapt();
  }

  return testGrid(grid, name);
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite test;

  // test with UG
  {
    using GridType = Dune::UGGrid<2>;
    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0}, {1.0,1.0}, {16,16});
    grid->setClosureType(GridType::ClosureType::NONE);
    test.subTest(testRefinement1(*grid, "UGGrid-1"));
  }
  {
    using GridType = Dune::UGGrid<2>;
    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0}, {2.0,1.0}, {2,1});
    grid->setClosureType(GridType::ClosureType::NONE);
    test.subTest(testRefinement2(*grid, "UGGrid-2"));
  }
  // test with ALU
#if HAVE_DUNE_ALUGRID
  {
    using GridType = Dune::ALUGrid<2,2, Dune::cube, Dune::nonconforming>;
    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid({0.0,0.0}, {1.0,1.0}, {16,16});
    test.subTest(testRefinement1(*grid, "ALUGrid-1"));
  }
  {
    using GridType = Dune::ALUGrid<2,2, Dune::cube, Dune::nonconforming>;
    auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0}, {2.0,1.0}, {2,1});
    grid->setClosureType(GridType::ClosureType::NONE);
    test.subTest(testRefinement2(*grid, "ALUGrid-2"));
  }
#endif

  return test.exit();
}
