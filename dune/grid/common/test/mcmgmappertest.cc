// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the MultipleCodimMultipleGeometryMapper
 */

#include <config.h>

#include <iostream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/uggrid.hh>
#include "../../../../doc/grids/gridfactory/hybridtestgrids.hh"

using namespace Dune;

/*!
 * \brief Layout template for edges and elements
 * This layout template is for use in the MultipleCodimMultipleGeomTypeMapper.
 * It selects edges and elements (entities with dim=1 or dim=dimgrid).
 *
 * \tparam dimgrid The dimension of the grid.
 */
template<int dimgrid>
struct MCMGElementEdgeLayout
{
  /*!
   * \brief: test whether entities of the given geometry type should be included in
   * the map
   */
  bool contains(GeometryType gt)
  {
    return (gt.dim() == dimgrid || gt.dim() == 1);
  }
};

/*!
 * \brief Check whether the index created for element data is unique,
 * consecutive and starting from zero.
 */
template <class Mapper, class GridView>
void checkElementDataMapper(const Mapper& mapper, const GridView& gridView)
{
  typedef typename GridView::template Codim<0>::Iterator Iterator;

  Iterator eIt    = gridView.template begin<0>();
  Iterator eEndIt = gridView.template end<0>();

  size_t min = 1;
  size_t max = 0;
  std::set<int> indices;

  for (; eIt!=eEndIt; ++eIt)
  {
    size_t index = mapper.index(*eIt);
    min = std::min(min, index);
    max = std::max(max, index);
    std::pair<std::set<int>::iterator, bool> status = indices.insert(index);

    if (!status.second)       // not inserted because already existing
      DUNE_THROW(GridError, "Mapper element index is not unique!");
  }

  if (min != 0)
    DUNE_THROW(GridError, "Mapper element index is not starting from zero!");

  if (max != gridView.indexSet().size(0) - 1)
    DUNE_THROW(GridError, "Mapper element index is not consecutive!");
}

/*!
 * \brief Check whether the index created for vertex data is consecutive
 * and starting from zero.
 */
template <class Mapper, class GridView>
void checkVertexDataMapper(const Mapper& mapper, const GridView& gridView)
{
  const size_t dim = GridView::dimension;
  typedef typename GridView::template Codim<0>::Iterator Iterator;

  Iterator eIt    = gridView.template begin<0>();
  Iterator eEndIt = gridView.template end<0>();

  size_t min = 1;
  size_t max = 0;
  std::set<int> indices;

  for (; eIt!=eEndIt; ++eIt)
  {
    size_t numVertices = eIt->subEntities(dim);
    for (size_t curVertex = 0; curVertex < numVertices; ++curVertex)
    {
      size_t index = mapper.subIndex(*eIt, curVertex, dim);
      min = std::min(min, index);
      max = std::max(max, index);
      indices.insert(index);
    }
  }

  if (min != 0)
    DUNE_THROW(GridError, "Mapper vertex index is not starting from zero!");

  if (max != gridView.indexSet().size(dim) - 1)
    DUNE_THROW(GridError, "Mapper vertex index is not consecutive!");

  for (size_t i = 0; i < max; ++i)
  {
    if (indices.find(i) == indices.end())
      DUNE_THROW(GridError, "Mapper vertex index is not consecutive!");
  }
}

/*!
 * \brief Check whether the index created for element and edge data is
 * unique, consecutive and starting from zero.
 */
template <class Mapper, class GridView>
void checkMixedDataMapper(const Mapper& mapper, const GridView& gridView)
{
  const size_t dim = GridView::dimension;
  typedef typename GridView::template Codim<0>::Iterator Iterator;

  Iterator eIt    = gridView.template begin<0>();
  Iterator eEndIt = gridView.template end<0>();

  size_t min = 1;
  size_t max = 0;
  std::set<int> indices;

  for (; eIt!=eEndIt; ++eIt)
  {
    // handle elements
    size_t index = mapper.index(*eIt);
    min = std::min(min, index);
    max = std::max(max, index);
    std::pair<std::set<int>::iterator, bool> status = indices.insert(index);

    if (!status.second)       // not inserted because already existing
      DUNE_THROW(GridError, "Mapper mixed index is not unique for elements!");

    // handle edges
    size_t numEdges = eIt->subEntities(dim-1);
    for (size_t curEdge = 0; curEdge < numEdges; ++curEdge)
    {
      index = mapper.subIndex(*eIt, curEdge, dim - 1);
      min = std::min(min, index);
      max = std::max(max, index);
      indices.insert(index);
    }
  }

  if (min != 0)
    DUNE_THROW(GridError, "Mapper mixed index is not starting from zero!");

  if (max != gridView.indexSet().size(0) + gridView.indexSet().size(dim - 1) - 1)
    DUNE_THROW(GridError, "Mapper mixed index is not consecutive!");

  for (size_t i = 0; i < max; ++i)
  {
    if (indices.find(i) == indices.end())
      DUNE_THROW(GridError, "Mapper mixed index is not consecutive!");
  }
}

/*!
 * \brief Run checks for a given grid.
 *
 * \param grid Grid to perform the checks.
 */
template<typename Grid>
void checkGrid(const Grid& grid)
{
  static const unsigned dim = Grid::dimension;
  // check element mapper
  // check leafMCMGMapper
  {   // check constructor without layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementLayout>
    leafMCMGMapper(grid);
    checkElementDataMapper(leafMCMGMapper, grid.leafGridView());
  }
  {   // check constructor with layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementLayout>
    leafMCMGMapper(grid, MCMGElementLayout<dim>());
    checkElementDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    {     // check constructor without layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementLayout>
      levelMCMGMapper(grid, i);
      checkElementDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementLayout>
      levelMCMGMapper(grid, i, MCMGElementLayout<dim>());
      checkElementDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
  }

  // check vertex mapper
  // check leafMCMGMapper
  {   // check constructor without layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGVertexLayout>
    leafMCMGMapper(grid);
    checkVertexDataMapper(leafMCMGMapper, grid.leafGridView());
  }
  {   // check constructor with layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGVertexLayout>
    leafMCMGMapper(grid, MCMGVertexLayout<dim>());
    checkVertexDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    {     // check constructor without layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGVertexLayout>
      levelMCMGMapper(grid, i);
      checkVertexDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGVertexLayout>
      levelMCMGMapper(grid, i, MCMGVertexLayout<dim>());
      checkVertexDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
  }

  // check mixed element and vertex mapper
  // check leafMCMGMapper
  {   // check constructor without layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementEdgeLayout>
    leafMCMGMapper(grid);
    checkMixedDataMapper(leafMCMGMapper, grid.leafGridView());
  }
  {   // check constructor with layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementEdgeLayout>
    leafMCMGMapper(grid, MCMGElementEdgeLayout<dim>());
    checkMixedDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    {     // check constructor without layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementEdgeLayout>
      levelMCMGMapper(grid, i);
      checkMixedDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid, MCMGElementEdgeLayout>
      levelMCMGMapper(grid, i, MCMGElementEdgeLayout<dim>());
      checkMixedDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
  }
}

int main(int argc, char** argv)
try
{
  // initialize MPI if necessary
  Dune::MPIHelper::instance(argc, argv);

  // Check grids with more than one element type.
  // So far only UGGrid does this, so we use it to test the mappers.

  //  Do the test for a 2d UGGrid
  {
    typedef UGGrid<2> Grid;
    std::unique_ptr<Grid> grid(make2DHybridTestGrid<Grid>());

    // create hybrid grid
    grid->mark(1, grid->leafGridView().begin<0>());
    grid->adapt();
    grid->globalRefine(1);

    checkGrid(*grid);
  }

  //  Do the  test for a 3d UGGrid
  {
    typedef UGGrid<3> Grid;
    std::unique_ptr<Grid> grid(make3DHybridTestGrid<Grid>());

    // create hybrid grid
    grid->mark(1, grid->leafGridView().begin<0>());
    grid->adapt();
    grid->globalRefine(1);

    checkGrid(*grid);
  }

  return EXIT_SUCCESS;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
