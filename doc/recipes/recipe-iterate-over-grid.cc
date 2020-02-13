/**
 * \page recipe-iterate-over-grid Iterating over a grid
 *
 * A fundamental task in many finite element-type discretizations is iterating
 * over the elements of a grid, for example in order to numerically compute integrals
 * on it. All DUNE grids provide a unified interface for this and related operations.
 *
 * First, we set up a simple structured grid.
 *
 * Then, we iterate over it
 * \snippet recipe-iterate-over-grid.cc iterate over grid view
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
  try{

    // Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // the simplest way to create a grid
    const int dim = 4;
    using Grid = Dune::YaspGrid<dim>;
    Dune::FieldVector<double,dim> len; for (auto& l:len) l=1.0;
    std::array<int,dim> cells; for (auto& c : cells) c=4;
    Grid grid(len,cells);

    // refine and extract grid view
    std::cout << "refining grid" << std::endl;
    grid.globalRefine(1);
    auto gv = grid.leafGridView();

    std::cout << "iterating over codim" << std::endl;
    const int codim = 2;
    for (const auto& e : entities(gv,Dune::Codim<codim>{}))
      if (!e.type().isCube()) std::cout << "not a cube" << std::endl;

    // [iterate over grid view]
    for (const auto& e : elements(gv)); // codim=0
    for (const auto& e : vertices(gv)); // codim=dim
    for (const auto& e : edges(gv));    // codim=dim-1
    for (const auto& e : facets(gv));   // codim=1
    //! [iterate over grid view]

    // access subentites
    std::cout << "iterating over subentities " << std::endl;
    for (const auto& e : elements(gv))
      for (unsigned int i=0; i<e.subEntities(codim); ++i)
        auto v = e.template subEntity<codim>(i);

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
