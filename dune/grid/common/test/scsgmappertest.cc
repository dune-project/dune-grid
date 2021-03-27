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
template <class Mapper, class GridView>
void checkElementDataMapper(const Mapper& mapper, const GridView& gridView)
{
  typedef typename GridView::IndexSet IndexSet;
  const IndexSet &is = gridView.indexSet();

  typedef typename Mapper::Index Index;

  typedef typename GridView::template Codim<0>::Iterator Iterator;

  Iterator eIt    = gridView.template begin<0>();
  Iterator eEndIt = gridView.template end<0>();

  Index min = 1;
  Index max = 0;
  std::set<Index> indices;
  typename Mapper::size_type size = mapper.size();

  if (size != is.size(0))
    DUNE_THROW(GridError, "Mapper size does not agree with index set "
               "size!");

  for (; eIt!=eEndIt; ++eIt) {

    Index index;
    bool contained = mapper.contains(*eIt, index);

    if ((is.types(0)[0] == eIt->type()) != contained)
      DUNE_THROW(GridError, "Mapper::contains() does not agree with the "
                 "element's geometry type!");
    if (!contained)
      continue;

    if (index != mapper.index(*eIt))
      DUNE_THROW(GridError, "Mapper::contains() and mapper.index() "
                 "compute different indices!");

    min = std::min(min, index);
    max = std::max(max, index);

    // typedef typename GridView::IntersectionIterator IntersectionIterator;

    // IntersectionIterator iIt    = gridView.ibegin(*eIt);
    // IntersectionIterator iEndIt = gridView.iend(*eIt);

    // for (; iIt!=iEndIt; ++iIt) {

    //     int oldindex = mapper.template map<1>(*eIt, iIt->numberInSelf());
    //     int index = mapper.map(*eIt, iIt->indexInInside(), 1);
    //     assert(oldindex == index);
    // }

    std::pair<typename std::set<Index>::iterator, bool> status = indices.insert(index);

    if (!status.second)       // not inserted because already existing
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

/*
   The MultipleGeometryMultipleCodimMapper only does something helpful on grids with more
   than one element type.  So far only UGGrids do this, so we use them to test the mapper.
 */

int main (int argc, char** argv) try
{
  // initialize MPI if necessary
  Dune :: MPIHelper::instance( argc, argv );

  // ////////////////////////////////////////////////////////////////////////
  //  Do the test for a 2d YaspGrid
  // ////////////////////////////////////////////////////////////////////////
  {
    static const int dim = 2;
    typedef YaspGrid<dim> GridType;
    typedef GridType::ctype ctype;

    Dune::FieldVector<ctype, dim> L(1.0);
    std::array<int, dim> s;
    std::fill(s.begin(), s.end(), 1);
    GridType grid(L,s);

    // create hybrid grid
    grid.mark(1, * grid.leafbegin<0>());
    grid.adapt();
    grid.globalRefine(1);

    LeafSingleCodimSingleGeomTypeMapper<GridType, 0> leafMapper(grid);
    std::vector<LevelSingleCodimSingleGeomTypeMapper<GridType, 0>> levelMappers;
    for (int i=0; i<=grid.maxLevel(); i++)
      levelMappers.push_back(LevelSingleCodimSingleGeomTypeMapper<GridType, 0>(grid, i));

    // Check mappers
    checkElementDataMapper(leafMapper, grid.leafGridView());
    for (std::size_t i=0; i<levelMappers.size(); i++)
      checkElementDataMapper(levelMappers[i], grid.levelGridView(i));

    // Refine grid and update mappers
    grid.globalRefine(1);
    leafMapper.update();
    for (std::size_t i=0; i<levelMappers.size(); i++)
      levelMappers[i].update();

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
