// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>
#include <numeric>

#include <dune/common/exceptions.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/globalindexset.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;

/** \brief Perform various consistency checks on a GlobalIndexSet object */
template <class GridView, int codim>
void checkIndexSet(const GridView& gridView,
                   const GlobalIndexSet<GridView>& indexSet)
{
  // Collect global indices on the current process
  std::vector<typename GlobalIndexSet<GridView>::Index> indices;
  for (auto it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it)
    // Loop over all subEntities
    for (size_t i=0; i<it->subEntities(codim); i++) {
      assert( indexSet.index(it->template subEntity<codim>(i)) == indexSet.subIndex(*it, i, codim) );
      indices.push_back(indexSet.index(it->template subEntity<codim>(i)));
    }

  /////////////////////////////////////////////////////////////////////////
  //  Sent the locally collected global indices to the root process
  /////////////////////////////////////////////////////////////////////////
  std::vector<int> sizes(gridView.comm().size());

  int share = indices.size();
  gridView.comm().template allgather<int>(&share, 1, sizes.data());

  std::vector<int> offsets(sizes.size());

  for (size_t k = 0; k < offsets.size(); ++k)
    offsets[k] = std::accumulate(sizes.begin(), sizes.begin() + k, 0);

  std::vector<typename GlobalIndexSet<GridView>::Index> indicesGlobal;
  if (gridView.comm().rank() == 0)
    indicesGlobal.resize(std::accumulate(sizes.begin(), sizes.end(), 0));

  gridView.comm().gatherv(indices.data(),
                          indices.size(),
                          indicesGlobal.data(),
                          sizes.data(),
                          offsets.data(),
                          0);   // root rank

  /////////////////////////////////////////////////////////////////////////////////
  //  Check whether the set of global indices is consecutive and starts at zero.
  //  (It may contain multiple entries, though.)
  /////////////////////////////////////////////////////////////////////////////////

  // To check we remove the duplicates
  std::sort(indicesGlobal.begin(), indicesGlobal.end());
  auto last = std::unique(indicesGlobal.begin(), indicesGlobal.end());
  indicesGlobal.erase(last, indicesGlobal.end());

  if (gridView.comm().rank()==0)
    for (size_t i=0; i<indicesGlobal.size(); i++)
      if ( (size_t) indicesGlobal[i] != i )
        DUNE_THROW(Exception, i << "th global index is not " << i);

}

int main(int argc, char* argv[]) try
{
#if HAVE_DUNE_UGGRID
  Dune::MPIHelper& mpiHelper =
#endif
  MPIHelper::instance(argc, argv);

  ////////////////////////////////////////////////////
  //  Create a distributed YaspGrid
  ////////////////////////////////////////////////////

#if HAVE_DUNE_UGGRID
  static const int dim = 2;
#endif
  // Disable the YaspGrid test for the time being.  It crashes in many situations
  // and I suspect that that's caused by bugs in YaspGrid.  I'll need to investigate that.
#if 0
  typedef YaspGrid<dim> GridType;

  std::array<int,dim> elements = {4, 4};
  FieldVector<double,dim> bbox = {10, 10};
  std::bitset<dim> periodic(0);
  unsigned int overlap = 1;

  GridType grid(MPI_COMM_WORLD, bbox, elements, periodic, overlap);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();
#endif

#if HAVE_DUNE_UGGRID
  typedef UGGrid<dim> GridType;

  std::array<unsigned int,dim> elements = { {8, 8} };
  FieldVector<double,dim> lower = {0, 0};
  FieldVector<double,dim> bbox = {10, 10};

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, bbox, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  grid->loadBalance();

  /////////////////////////////////////////////////////
  //  Create and check global index sets
  /////////////////////////////////////////////////////

  // elements
  if (mpiHelper.rank() == 0)
    std::cout << "Elements" << std::endl;
  GlobalIndexSet<GridView> elementIndexSet(gridView,0);
  checkIndexSet<GridView,0>(gridView, elementIndexSet);

  // edges
  if (mpiHelper.rank() == 0)
    std::cout << "Edges" << std::endl;
  GlobalIndexSet<GridView> edgeIndexSet(gridView,1);
  checkIndexSet<GridView,1>(gridView, edgeIndexSet);

  // vertices
  if (mpiHelper.rank() == 0)
    std::cout << "Vertices" << std::endl;
  GlobalIndexSet<GridView> vertexIndexSet(gridView,2);
  checkIndexSet<GridView,2>(gridView, vertexIndexSet);
#endif

  return 0;

}
catch (const Dune::Exception &e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
