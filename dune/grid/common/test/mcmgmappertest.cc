// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#include "../../../../doc/grids/gridfactory/hybridtestgrids.hh"
#endif

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
template <class G, class M, class I, class GridView>
void checkElementDataMapper(const Dune::Mapper<G, M, I>& mapper, const GridView& gridView)
{
  using Mapper = Dune::Mapper<G, M, I>;
  using Index = typename Mapper::Index;
  Index min = 1;
  Index max = 0;
  std::set<int> indices;

  for (const auto& element : elements(gridView))
  {
    const auto index = mapper.index(element);
    min = std::min(min, index);
    max = std::max(max, index);
    [[maybe_unused]] const auto [it, wasInserted] = indices.insert(index);

    if (!wasInserted) // not inserted because already existing
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
template <class G, class M, class I, class GridView>
void checkVertexDataMapper(const Dune::Mapper<G, M, I>& mapper, const GridView& gridView)
{
  const size_t dim = GridView::dimension;
  using Mapper = Dune::Mapper<G, M, I>;
  using Index = typename Mapper::Index;

  Index min = 1;
  Index max = 0;
  std::set< Index > indices;

  for (const auto& element : elements(gridView))
  {
    size_t numVertices = element.subEntities(dim);
    for (size_t curVertex = 0; curVertex < numVertices; ++curVertex)
    {
      Index testIdx = Index(-1);
      bool contains = mapper.contains(element, curVertex, dim, testIdx );
      Index index = mapper.subIndex(element, curVertex, dim);

      if( contains && (index != testIdx) )
      {
        DUNE_THROW(GridError, "subIndex and contains do not return the same index!");
      }

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
 * \note Tests the multiple indices per entity interface (types/indices) specific to the MCMGMapper
 */
template <class Mapper, class GridView>
void checkMixedDataMapper(const Mapper& mapper, const GridView& gridView)
{
  using Index = typename Mapper::Index;
  const size_t dim = GridView::dimension;

  Index min = 1;
  Index max = 0;
  std::set<int> indices;

  Index elementBlockSize = 0, edgeBlockSize = 0;

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
      [[maybe_unused]] const auto [it, wasInserted] = indices.insert(i);

      if (!wasInserted) // not inserted because already existing
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

      Index testIdx = Index(-1);
      const bool contains = mapper.contains(element, curEdge, dim-1, testIdx );
      const Index index = mapper.subIndex(element, curEdge, dim-1);

      if( contains && (index != testIdx) )
      {
        DUNE_THROW(GridError, "subIndex and contains do not return the same index!");
      }
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

template <class G, class M, class I, class GV>
void update(Dune::Mapper<G, M, I>& mapper, const GV& gv)
{ mapper.update(gv); }

/*!
 * \brief Run checks for a given grid.
 *
 * \param grid Grid to perform the checks.
 */
template<typename Grid>
void checkGrid(Grid& grid)
{
  using LeafMCMGMapper = MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView>;
  using LevelMCMGMapper = MultipleCodimMultipleGeomTypeMapper<typename Grid::LevelGridView>;

  // Create various mappers *******************************************************************

  // Create element mappers for leaf and levels
  LeafMCMGMapper leafElementMCMGMapper(grid.leafGridView(), mcmgElementLayout());
  std::vector<LevelMCMGMapper> levelElementMCMGMappers;
  for (int i = 0; i <= grid.maxLevel(); i++)
    levelElementMCMGMappers.emplace_back(grid.levelGridView(i), mcmgElementLayout());

  // Create vertex mappers for leaf and levels
  LeafMCMGMapper leafVertexMCMGMapper(grid.leafGridView(), mcmgVertexLayout());
  std::vector<LevelMCMGMapper> levelVertexMCMGMappers;
  for (int i = 0; i <= grid.maxLevel(); i++)
    levelVertexMCMGMappers.emplace_back(grid.levelGridView(i), mcmgVertexLayout());

  // Create mixed element and edge mappers for leaf and levels
  const auto elementEdgeLayout = [](GeometryType gt, unsigned int dimgrid) {
    return (gt.dim() == dimgrid)? 3 : (gt.dim() == 1? 2 : 0);
  };
  LeafMCMGMapper leafMixedMCMGMapper(grid.leafGridView(), elementEdgeLayout);
  std::vector<LevelMCMGMapper> levelMixedMCMGMappers;
  for (int i = 0; i <= grid.maxLevel(); i++)
    levelMixedMCMGMappers.emplace_back(grid.levelGridView(i), elementEdgeLayout);

  // Check mappers *******************************************************************

  checkElementDataMapper(leafElementMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelElementMCMGMappers.size(); i++)
    checkElementDataMapper(levelElementMCMGMappers[i], grid.levelGridView(i));

  checkVertexDataMapper(leafVertexMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelVertexMCMGMappers.size(); i++)
    checkVertexDataMapper(levelVertexMCMGMappers[i], grid.levelGridView(i));

  checkMixedDataMapper(leafMixedMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelMixedMCMGMappers.size(); i++)
    checkMixedDataMapper(levelMixedMCMGMappers[i], grid.levelGridView(i));

  // Refine grid and update mappers *******************************************************************

  grid.globalRefine(1);

  update(leafElementMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelElementMCMGMappers.size(); i++)
    update(levelElementMCMGMappers[i], grid.levelGridView(i));

  update(leafVertexMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelVertexMCMGMappers.size(); i++)
    update(levelVertexMCMGMappers[i], grid.levelGridView(i));

  update(leafMixedMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelMixedMCMGMappers.size(); i++)
    update(levelMixedMCMGMappers[i], grid.levelGridView(i));

  // Check mappers *******************************************************************

  checkElementDataMapper(leafElementMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelElementMCMGMappers.size(); i++)
    checkElementDataMapper(levelElementMCMGMappers[i], grid.levelGridView(i));

  checkVertexDataMapper(leafVertexMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelVertexMCMGMappers.size(); i++)
    checkVertexDataMapper(levelVertexMCMGMappers[i], grid.levelGridView(i));

  checkMixedDataMapper(leafMixedMCMGMapper, grid.leafGridView());
  for (std::size_t i = 0; i < levelMixedMCMGMappers.size(); i++)
    checkMixedDataMapper(levelMixedMCMGMappers[i], grid.levelGridView(i));

}

int main(int argc, char** argv)
try
{
  // initialize MPI if necessary
  Dune::MPIHelper::instance(argc, argv);

  // Check grids with more than one element type.
  // So far only UGGrid does this, so we use it to test the mappers.

#if HAVE_DUNE_UGGRID
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
#endif

  //  Do the test for a 2d YaspGrid
  {
    Dune::FieldVector< double, 2 > lower( 0 );
    Dune::FieldVector< double, 2 > upper( 1 );
    std::array<unsigned int,2 > elements = {{ 4 ,4 }};

    typedef Dune::YaspGrid<2> Grid;
    auto grid = Dune::StructuredGridFactory< Grid >::createCubeGrid( lower, upper, elements );
    // create hybrid grid
    checkGrid(*grid);
  }

  //  Do the test for a 3d YaspGrid
  {
    Dune::FieldVector< double, 3 > lower( 0 );
    Dune::FieldVector< double, 3 > upper( 1 );
    std::array<unsigned int,3 > elements = {{ 4 ,4, 4 }};

    typedef Dune::YaspGrid<3> Grid;
    auto grid = Dune::StructuredGridFactory< Grid >::createCubeGrid( lower, upper, elements );
    // create hybrid grid
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
