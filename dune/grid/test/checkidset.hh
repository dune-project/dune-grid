// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKIDSET_HH
#define DUNE_GRID_TEST_CHECKIDSET_HH

#include <map>
#include <utility>

#include <dune/common/float_cmp.hh>
#include <dune/common/hybridutilities.hh>

/** @file
   @brief Unit tests for IdSet implementations
 */


namespace Dune
{

  template <class Grid, class IdSet>
  void checkInjectivity(const Grid& grid, const IdSet& idSet)
  {
    const int dim = Grid::dimension;

    using IdType = typename IdSet::IdType;
    using GlobalCoordinate = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;

    std::map<IdType,GlobalCoordinate> idContainer;

    // Loop over all entities of all codimensions on all grid levels
    for (int i=0; i<=grid.maxLevel(); ++i)
    {
      const auto gridView = grid.levelGridView(i);

      for (const auto& element: elements(gridView))
      {
        // Loop over all faces of all codimensions
        Hybrid::forEach( std::make_index_sequence< dim+1 >{},
                         [ & ]( auto codim )
        {
          for (size_t face=0; face<element.subEntities(codim); ++face)
          {
            const auto& entity = element.template subEntity<codim>(face);

            auto id = idSet.id(entity);

            // Has the same id already been used by a different entity?
            if (idContainer.find(id) != idContainer.end())
            {
              // Yes.  Then either we have seen the same entity before, or we are now
              // on the copy of an entity we have seen before.  In either case we must
              // have the same entity center.
              // CAVEAT: This last reasoning does not hold if the grid uses parametrized
              // elements or parametrized boundaries.
              if (! FloatCmp::eq(entity.geometry().center(), idContainer[id]))
                DUNE_THROW(GridError, "IdSet is not injective");
            }
            else
            {
              idContainer[id] = entity.geometry().center();
            }
          }
        });
      }
    }
  }


  // Run all available tests for a given IdSet
  template <class Grid, class IdSet>
  void checkIdSet ( const Grid &grid, const IdSet& idSet)
  {
    checkInjectivity(grid, idSet);
  }

} // end namespace Dune

#endif // DUNE_GRID_TEST_CHECKIDSET_HH
