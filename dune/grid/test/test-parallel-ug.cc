// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Test parallel interface if a parallel UG is used

#include <config.h>

#include <dune/grid/test/checkparallelug.hh>

template <int dim>
void test(bool localRefinement) {

  typedef Dune::UGGrid<dim> GridType;

  // Create grid
  Dune::StructuredGridFactory<GridType> structuredGridFactory;

  Dune::FieldVector<double,dim> lowerLeft(0);
  Dune::FieldVector<double,dim> upperRight(1);
  std::array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 4);
  std::shared_ptr<GridType> grid = structuredGridFactory.createCubeGrid(lowerLeft, upperRight, elements);

  // Balance (with check) and check the resulting grid
  Dune::GridCheck::testParallelUGLoadBalance (grid);
  Dune::GridCheck::testParallelUG(grid);

  // Refinement
  if (!localRefinement)
    grid->globalRefine(1);
  else {
    // mark all elements with x-coordinate < 0.5 for refinement
    typedef typename GridType::LeafGridView LeafGV;
    typename LeafGV::template Codim<0>::Iterator
    it = grid->leafGridView().template begin<0>();
    const typename LeafGV::template Codim<0>::Iterator
    &endIt = grid->leafGridView().template end<0>();
    for (; it != endIt; ++it) {
      int nRefine = 1;
      if (it->geometry().center()[0] < 0.5)
        grid->mark(nRefine, *it);
    }

    // adapt the grid
    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();
  }

  // And check again
  Dune::GridCheck::testParallelUG(grid);

}

int main (int argc , char **argv) try
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper &mpiHelper =
    Dune::MPIHelper::instance(argc, argv);

  std::cout << "This is process "
            << mpiHelper.rank() + 1
            << " of "
            << mpiHelper.size()
            << ", PID "
            << getpid()
            << " .\n";

  // test 2D/3D grids with uniform/local refinement
  test<2>(false);
  test<2>(true);
  test<3>(false);
  test<3>(true);

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
