// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the SingleCodimSingleGeomTypeMapper
 */

#include <config.h>

#include <iostream>
#include <set>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

using namespace Dune;

// /////////////////////////////////////////////////////////////////////////////////
//   Check whether the index created for element data is unique, consecutive
//   and starting from zero.
// /////////////////////////////////////////////////////////////////////////////////
template <class G, class M, class I, class GridView>
void checkElementDataMapper(const Dune::Mapper<G, M, I>& mapper, const GridView& gridView)
{
  const auto& is = gridView.indexSet();

  using Mapper = Dune::Mapper<G, M, I>;
  using Index = typename Mapper::Index;

  Index min = 1;
  Index max = 0;
  std::set<Index> indices;
  const auto size = mapper.size();

  if (size != is.size(0))
    DUNE_THROW(GridError, "Mapper size does not agree with index set "
               "size!");

  for (const auto& element : elements(gridView))
  {
    Index index;
    bool contained = mapper.contains(element, index);

    if ((is.types(0)[0] == element.type()) != contained)
      DUNE_THROW(GridError, "Mapper::contains() does not agree with the "
                 "element's geometry type!");
    if (!contained)
      continue;

    if (index != mapper.index(element))
      DUNE_THROW(GridError, "Mapper::contains() and mapper.index() "
                 "compute different indices!");

    min = std::min(min, index);
    max = std::max(max, index);

    // typedef typename GridView::IntersectionIterator IntersectionIterator;

    // IntersectionIterator iIt    = gridView.ibegin(element);
    // IntersectionIterator iEndIt = gridView.iend(element);

    // for (; iIt!=iEndIt; ++iIt) {

    //     int oldindex = mapper.template map<1>(element, iIt->numberInSelf());
    //     int index = mapper.map(element, iIt->indexInInside(), 1);
    //     assert(oldindex == index);
    // }

    [[maybe_unused]] const auto [it, wasInserted] = indices.insert(index);

    if (!wasInserted)       // not inserted because already existing
      DUNE_THROW(GridError, "Mapper element index is not unique!");
  }

  if (size!=indices.size())
    DUNE_THROW(GridError, "Mapper size does not agree with the number of "
               "elements counted by iteration.");

  if (min!=0)
    DUNE_THROW(GridError, "Mapper element index is not starting from zero!");

  if (max+1!=size)
    DUNE_THROW(GridError, "Mapper element index is not consecutive!");

}

template <class G, class M, class I, class GV>
void update(Dune::Mapper<G, M, I>& mapper, const GV& gv)
{ mapper.update(gv); }

int main (int argc, char** argv) try
{
  // initialize MPI if necessary
  Dune :: MPIHelper::instance( argc, argv );

  // ////////////////////////////////////////////////////////////////////////
  //  Do the test for a 2d YaspGrid
  // ////////////////////////////////////////////////////////////////////////
  {
    static const int dim = 2;
    using Grid = YaspGrid<dim>;
    using ctype = Grid::ctype;

    Dune::FieldVector<ctype, dim> L(1.0);
    std::array<int, dim> s;
    std::fill(s.begin(), s.end(), 1);
    Grid grid(L, s);

    // create hybrid grid
    grid.mark(1, * grid.leafbegin<0>());
    grid.adapt();
    grid.globalRefine(1);

    using LeafGridView = typename Grid::LeafGridView;
    using LeafMapper = SingleCodimSingleGeomTypeMapper<LeafGridView, 0>;
    LeafMapper leafMapper(grid.leafGridView());

    using LevelGridView = typename Grid::LevelGridView;
    using LevelMapper = SingleCodimSingleGeomTypeMapper<LevelGridView, 0>;
    std::vector<LevelMapper> levelMappers;
    for (int i=0; i<=grid.maxLevel(); i++)
      levelMappers.emplace_back(grid.levelGridView(i));

    // Check mappers
    checkElementDataMapper(leafMapper, grid.leafGridView());
    for (std::size_t i=0; i<levelMappers.size(); i++)
      checkElementDataMapper(levelMappers[i], grid.levelGridView(i));

    // Refine grid and update mappers
    grid.globalRefine(1);
    update(leafMapper, grid.leafGridView());
    for (std::size_t i=0; i<levelMappers.size(); i++)
      update(levelMappers[i], grid.levelGridView(i));

    // Check mappers
    checkElementDataMapper(leafMapper, grid.leafGridView());
    for (std::size_t i=0; i<levelMappers.size(); i++)
      checkElementDataMapper(levelMappers[i], grid.levelGridView(i));

  }

  return 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
