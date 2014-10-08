// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/globalindexset.hh>

using namespace Dune;

/** \brief Perform various consistency checks on a GlobalIndexSet object */
template <class GridView, int codim>
void checkIndexSet(const GridView& gridView,
                   const GlobalIndexSet<GridView,codim>& indexSet)
{
  // Collect global indices on the current process
  std::vector<typename GlobalIndexSet<GridView,codim>::Index> indices;
  for (auto it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it)
    // Loop over all subEntities
    for (size_t i=0; i<it->subEntities(codim); i++)
      indices.push_back(indexSet.globalIndex(*it->template subEntity<codim>(i)));

  /////////////////////////////////////////////////////////////////////////
  //  Sent the locally collected global indices to the root process
  /////////////////////////////////////////////////////////////////////////
  std::vector<int> sizes(gridView.comm().size());

  int share = indices.size();
  gridView.comm().template allgather<int>(&share, 1, sizes.data());

  std::vector<int> offsets(sizes.size());

  for (size_t k = 0; k < offsets.size(); ++k)
    offsets[k] = std::accumulate(sizes.begin(), sizes.begin() + k, 0);

  std::vector<typename GlobalIndexSet<GridView,codim>::Index> indicesGlobal;
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
  std::sort(indicesGlobal.begin(), indicesGlobal.end());

  if (gridView.comm().rank()==0)
  {
    if (indicesGlobal[0] != 0)
      DUNE_THROW(Exception, "Global index set does not contain the index '0'");

    for (size_t i=0; i<indicesGlobal.size()-1; i++)
      if (indicesGlobal[i+1] > indicesGlobal[i]+1)
        DUNE_THROW(Exception, "Global index set does not contain the index '" << indicesGlobal[i]+1 << "'");
  }
}

int main(int argc, char* argv[]) try
{
  Dune::MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);

  ////////////////////////////////////////////////////
  //  Create a distributed YaspGrid
  ////////////////////////////////////////////////////

  static const int dim = 2;
  typedef YaspGrid<dim> GridType;

  array<int,dim> elements = {4, 4};
  FieldVector<double,dim> bbox = {10, 10};
  std::bitset<dim> periodic(0);
  uint overlap = 1;

  GridType grid(MPI_COMM_WORLD, bbox, elements, periodic, overlap);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////
  //  Create and check global index sets
  /////////////////////////////////////////////////////

  // elements
  GlobalIndexSet<GridView,0> elementIndexSet(gridView);
  checkIndexSet(gridView, elementIndexSet);

  // edges
  GlobalIndexSet<GridView,1> edgeIndexSet(gridView);
  checkIndexSet(gridView, edgeIndexSet);

  // vertices
  GlobalIndexSet<GridView,2> vertexIndexSet(gridView);
  checkIndexSet(gridView, vertexIndexSet);

  return 0;

}
catch (const Dune::Exception &e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
