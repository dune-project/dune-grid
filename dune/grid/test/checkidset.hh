// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKIDSET_HH
#define DUNE_GRID_TEST_CHECKIDSET_HH

#include <iostream>
#include <map>
#include <utility>

#include <dune/common/float_cmp.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/exceptions.hh>

/** @file
   @brief Unit tests for IdSet implementations
 */


namespace Dune
{
  // Check the type used for ids
  template <class Grid, class IdSet>
  void checkIdType(const Grid& grid, const IdSet& idSet)
  {
    auto gridView = grid.leafGridView();
    auto begin = gridView.template begin<0>();

    // skip test for empty grids
    if( begin == gridView.template end<0>() )
      return ;

    // Get some entity for testing
    auto entity = *begin;

    // The IdSet class exports a type IdSet, and this is the type actually
    // used for return values of the id class
    static_assert(std::is_same<typename IdSet::IdType, decltype(idSet.id(entity))>::value,
                  "IdSet::IdType does not match the return value of idSet.id()");
    static_assert(std::is_same<typename IdSet::IdType, decltype(idSet.subId(entity,0,0))>::value,
                  "IdSet::IdType does not match the return value of idSet.subId()");

    // Get some id for testing
    auto id = idSet.id(entity);

    // Test for operator<
    if (id<id)
      DUNE_THROW(GridError, "operator< does not implement an ordering");

    // Test for operator==
    if (!(id==id))
      DUNE_THROW(GridError, "operator== is not symmetric");

    // Test for operator!=
    if (id!=id)
      DUNE_THROW(GridError, "operator!= is not the negation of operator==");

    // Test for the default constructor
    typename IdSet::IdType someId;

    // Test for the copy constructor
    typename IdSet::IdType someOtherId(idSet.id(entity));

    // copy assignment
    someId = someOtherId;

    // Output operator
    std::cout << someId << std::endl;
  }

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
              if (! FloatCmp::eq(entity.geometry().center(), idContainer[id], 1e-12 ))
                DUNE_THROW(GridError, "IdSet is not injective");
            }
            else
            {
              idContainer[id] = entity.geometry().center();
            }

            // While we are here: Do subEntity.id and subId return the same value?
            if (id != idSet.subId(element,face,codim))
              DUNE_THROW(GridError, "subEntity.id and subId do not return the same value!");
          }
        });
      }
    }
  }


  // Run all available tests for a given IdSet
  template <class Grid, class IdSet>
  void checkIdSet ( const Grid &grid, const IdSet& idSet)
  {
    checkIdType(grid,idSet);
    checkInjectivity(grid, idSet);
  }

} // end namespace Dune

#endif // DUNE_GRID_TEST_CHECKIDSET_HH
