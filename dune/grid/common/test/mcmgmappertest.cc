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
  size_t min = 1;
  size_t max = 0;
  std::set<int> indices;

  for (const auto& element : elements(gridView))
  {
    size_t index = mapper.index(element);
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

  size_t min = 1;
  size_t max = 0;
  std::set<int> indices;

  for (const auto& element : elements(gridView))
  {
    size_t numVertices = element.subEntities(dim);
    for (size_t curVertex = 0; curVertex < numVertices; ++curVertex)
    {
      size_t index = mapper.subIndex(element, curVertex, dim);
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
  typedef typename Mapper::Index Index;
  const size_t dim = GridView::dimension;

  Index min = 1;
  Index max = 0;
  std::set<int> indices;

  Index elementBlockSize = 0, edgeBlockSize = 0;
  Index index;

  for (const auto& element : elements(gridView))
  {
    const auto &gts = mapper.types(0);
    if (std::find(gts.begin(),gts.end(),element.type()) == gts.end())
      DUNE_THROW(GridError, "Mapper mixed index does not have correct geometry type information");

    auto block = mapper.indices(element);
    if (block.empty())
      DUNE_THROW(GridError, "Mapper mixed index does not have element indices");
    if (elementBlockSize == 0) elementBlockSize = block.size();
    if (elementBlockSize != block.size())
      DUNE_THROW(GridError, "Mapper mixed index does not have the same block size on all elements");
    if ( block.size() != mapper.size(element.type()) )
      DUNE_THROW(GridError, "Mapper mixed index does not return correct size for given geometry type ");

    // handle elements
    min = std::min(min, *(block.begin()));
    max = std::max(max, *(block.end())-1);
    for (Index i : block)
    {
      std::pair<std::set<int>::iterator, bool> status = indices.insert(i);

      if (!status.second)       // not inserted because already existing
        DUNE_THROW(GridError, "Mapper mixed index is not unique for elements!");
    }

    // handle edges
    size_t numEdges = element.subEntities(dim-1);
    for (size_t curEdge = 0; curEdge < numEdges; ++curEdge)
    {
      auto block = mapper.indices(element,curEdge,dim-1);
      if (block.empty())
        DUNE_THROW(GridError, "Mapper mixed index does not have edges indices");
      if (edgeBlockSize == 0) edgeBlockSize = block.size();
      else if (edgeBlockSize != block.size())
        DUNE_THROW(GridError, "Mapper mixed index does not have the same block size on all edges");
      min = std::min(min, *(block.begin()));
      max = std::max(max, *(block.end())-1);
      for (Index i : block)
        indices.insert(i);
    }
  }

  if (min != 0)
    DUNE_THROW(GridError, "Mapper mixed index is not starting from zero!");

  if (max != gridView.indexSet().size(0)*elementBlockSize +
             gridView.indexSet().size(dim - 1)*edgeBlockSize - 1)
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
  DUNE_NO_DEPRECATED_BEGIN
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
  DUNE_NO_DEPRECATED_END
  {
    LeafMultipleCodimMultipleGeomTypeMapper<Grid>
    leafMCMGMapper(grid, mcmgElementLayout());
    checkElementDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    DUNE_NO_DEPRECATED_BEGIN
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
    DUNE_NO_DEPRECATED_END
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid>
      levelMCMGMapper(grid, i, mcmgElementLayout());
      checkElementDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
  }

  // check vertex mapper
  // check leafMCMGMapper
  DUNE_NO_DEPRECATED_BEGIN
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
  DUNE_NO_DEPRECATED_END
  {
    LeafMultipleCodimMultipleGeomTypeMapper<Grid>
    leafMCMGMapper(grid, mcmgVertexLayout());
    checkVertexDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    DUNE_NO_DEPRECATED_BEGIN
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
    DUNE_NO_DEPRECATED_END
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid>
      levelMCMGMapper(grid, i, mcmgVertexLayout());
      checkVertexDataMapper(levelMCMGMapper, grid.levelGridView(i));
    }
  }

  // check mixed element and edge mapper
  const auto elementEdgeLayout = [](GeometryType gt, int dimgrid) {
    return (gt.dim() == dimgrid)? 3 : (gt.dim() == 1? 2 : 0);
  };
  // check leafMCMGMapper
  DUNE_NO_DEPRECATED_BEGIN
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
  DUNE_NO_DEPRECATED_END
  {   // check constructor with layout class
    LeafMultipleCodimMultipleGeomTypeMapper<Grid>
    leafMCMGMapper(grid, elementEdgeLayout);
    checkMixedDataMapper(leafMCMGMapper, grid.leafGridView());
  }

  // check levelMCMGMapper
  for (int i = 2; i <= grid.maxLevel(); i++)
  {
    DUNE_NO_DEPRECATED_BEGIN
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
    DUNE_NO_DEPRECATED_END
    {     // check constructor with layout class
      LevelMultipleCodimMultipleGeomTypeMapper<Grid>
      levelMCMGMapper(grid, i, elementEdgeLayout);
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
    grid->mark(1, *grid->leafGridView().begin<0>());
    grid->adapt();
    grid->globalRefine(1);

    checkGrid(*grid);
  }

  //  Do the  test for a 3d UGGrid
  {
    typedef UGGrid<3> Grid;
    std::unique_ptr<Grid> grid(make3DHybridTestGrid<Grid>());

    // create hybrid grid
    grid->mark(1, *grid->leafGridView().begin<0>());
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
