// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file COPYING in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
/**
 * \page recipe-iterate-over-grid Iterating over a grid
 *
 * A fundamental task in many finite element-type discretizations is iterating
 * over the elements of a grid, for example in order to numerically compute integrals
 * on it. All DUNE grids provide a unified interface for this and related operations.
 *
 * First, we set up a simple structured grid. Here we use the \ref Dune::YaspGrid class in
 * four dimensions in order to demonstrate that even dimension larger than three
 * can be used.
 * \snippet recipe-iterate-over-grid.cc set up grid
 *
 * Grids in Dune are hierarchical, i.e. they are organized into levels
 * (originating from refinement) and entities that are not further refined.
 * Each of these subsets is accessible via Dune::GridView and iteration over
 * is possible only over grid views. So we extract the Dune::LeafGridView:
 * \snippet recipe-iterate-over-grid.cc extract gridview
 *
 * Now we can iterate over all the entities of a certain codimension in a grid view
 * using a range-based for loop:
 * \snippet recipe-iterate-over-grid.cc iterate over codim
 * As an example we extract the type of an element and test for it to be a cube.
 *
 * Instead of specifying the codimension explicitly you can also
 * use the following predefined names in the range-based for loop:
 * \snippet recipe-iterate-over-grid.cc iterate over grid view
 *
 * Finally, from entities of codimension 0 (aka elements) you can access
 * all subentities of all codimensions using the following code:
 * \snippet recipe-iterate-over-grid.cc access to subentities
 *
 * Full example code: @ref recipe-iterate-over-grid.cc
 * \example recipe-iterate-over-grid.cc
 * See explanation at @ref recipe-iterate-over-grid
 */


// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-grid includes
#include <dune/grid/yaspgrid.hh>


int main(int argc, char** argv)
{
  // Maybe initialize Mpi
  [[maybe_unused]] Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // [set up grid]
  const int dim = 4;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> len; for (auto& l : len) l=1.0;
  std::array<int,dim> cells; for (auto& c : cells) c=4;
  Grid grid(len,cells);
  //! [set up grid]

  // [extract gridview]
  auto gv = grid.leafGridView();
  //! [extract gridview]

  // [iterate over codim]
  const int codim = 2;
  for (const auto& e : entities(gv,Dune::Codim<codim>{}))
    if (!e.type().isCube()) std::cout << "not a cube" << std::endl;
  //! [iterate over codim]

  // [iterate over grid view]
  for ([[maybe_unused]] const auto& e : elements(gv))
    ; // codim=0
  for ([[maybe_unused]] const auto& e : vertices(gv))
    ; // codim=dim
  for ([[maybe_unused]] const auto& e : edges(gv))
    ;    // codim=dim-1
  for ([[maybe_unused]] const auto& e : facets(gv))
    ;   // codim=1
  //! [iterate over grid view]

  // [access to subentities]
  const int mycodim = 2;
  for (const auto& e : elements(gv))
    for (unsigned int i=0; i<e.subEntities(mycodim); ++i)
      [[maybe_unused]] auto v = e.template subEntity<codim>(i);
  //! [access to subentities]
}
