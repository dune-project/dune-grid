// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKADAPTATION_HH
#define DUNE_GRID_TEST_CHECKADAPTATION_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

#include <dune/grid/common/capabilities.hh>

/** \file
    \brief A test for the Adaptation interface
 */

template <class GridType>
void markForAdaptation(GridType& grid, const int marker )
{
  using namespace Dune;

  // Loop over all levels except the lowest one
  for (int level = 0; level <= grid.maxLevel(); ++level)
  {
    typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
    typedef typename ElementIterator :: Entity EntityType ;

    const ElementIterator eEndIt = grid.levelGridView(level).template end<0>();
    for (ElementIterator it = grid.levelGridView(level).template begin<0>();
         it != eEndIt; ++it )
    {
      const EntityType& entity = *it ;

      // marking should only be possible on leaf entities
      const bool marked = grid.mark(marker, entity );

      if( marker < 0 && entity.level() == 0 )
      {
        if( marked )
        {
          DUNE_THROW(InvalidStateException,"Macro entity was marked for coarsening");
        }
      }
      else
      {
        if( ! entity.isLeaf() && marked == true )
        {
          DUNE_THROW(InvalidStateException,"Non-leaf entity was marked");
        }

        if( entity.isLeaf() && marked == false )
        {
          DUNE_THROW(InvalidStateException,"Leaf entity could not be marked");
        }
      }
    }
  }
}

template <class EntityType>
void checkHierarchy(const EntityType& entity)
{
  using namespace Dune;

  if( ! entity.isLeaf() )
  {
    // only leaf entities can vanish
    if( entity.mightVanish() )
    {
      DUNE_THROW(InvalidStateException,"Non-leaf entity marked might vanish");
    }

    const int level = entity.level() + 1 ;

    typedef typename EntityType :: HierarchicIterator HierarchicIterator;
    const HierarchicIterator end = entity.hend( level );
    for(HierarchicIterator it = entity.hbegin( level ); it != end; ++it)
    {
      checkHierarchy( *it );
    }
  }
  else
  {
    if( entity.level() == 0 )
    {
      if( entity.mightVanish() )
        DUNE_THROW(InvalidStateException,"entity.mightVanish() returns true for macro element");
    }
    else
    {
      // on level 0 this should be false
      if( ! entity.mightVanish() )
      {
        DUNE_THROW(InvalidStateException,"entity.mightVanish() returns wrong result");
      }
    }
  }
}

template <class GridType>
void checkAdaptRefinement(GridType& grid, const bool greenClosure = false )
{
  using namespace Dune;

  // skip empty grids
  if ( grid.levelGridView(0).template begin<0>() == grid.levelGridView(0).template end<0>() )
    return;

  // some things are different for bisection grids
  const bool bisectionGrid =
    Capabilities::isLevelwiseConforming<GridType> :: v  == false &&
    Capabilities::isLeafwiseConforming<GridType>  :: v  == true ;

  for(int i=0; i<2; ++i)
  {
    const int oldMaxLevel = grid.maxLevel();

    // mark all leaf entities for refinement
    markForAdaptation( grid , 1 );

    bool markedCoarsen = grid.preAdapt();
    if( markedCoarsen != greenClosure )
    {
      DUNE_THROW(InvalidStateException,"grid.preAdapt() does not return correct information");
    }

    bool refined = grid.adapt() ;
    if( ! refined )
    {
      DUNE_THROW(InvalidStateException,"grid.adapt() returns wrong information");
    }

    /// check new max level
    if( grid.maxLevel () <= oldMaxLevel )
    {
      DUNE_THROW(InvalidStateException,"grid.maxLevel() wrong after refinement " << oldMaxLevel << " " << grid.maxLevel() );
    }

    // Loop over all levels except the lowest one
    for (int level = 0 ; level <= grid.maxLevel(); ++level )
    {
      typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
      typedef typename ElementIterator :: Entity EntityType ;
      ElementIterator eEndIt = grid.levelGridView(level).template end<0>();

      for (ElementIterator it = grid.levelGridView(level).template begin<0>();
           it != eEndIt; ++ it)
      {
        const EntityType& entity = *it ;

        // this check fails on biscetion grids
        if( ! bisectionGrid )
        {
          if( entity.isLeaf () != entity.isNew () )
          {
            DUNE_THROW(InvalidStateException,"isNew information on entity gives wrong result");
          }
        }
        else  // at least all leafs have to be new
        {
          if( entity.isLeaf () && ! entity.isNew () )
          {
            DUNE_THROW(InvalidStateException,"isNew information on entity gives wrong result");
          }
        }
      }
    }

    grid.postAdapt();

    // Loop over all levels except the lowest one
    for (int level = 0 ; level <= grid.maxLevel(); ++level )
    {
      typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
      ElementIterator eEndIt = grid.levelGridView(level).template end<0>();

      for (ElementIterator it = grid.levelGridView(level).template begin<0>();
           it != eEndIt; ++ it)
      {
        if( it->isNew () )
        {
          DUNE_THROW(InvalidStateException,"After postAdapt() was called no entity is new, i.e., isNew() == false");
        }
      }
    }
  }
}

template <class GridType>
void checkAdaptation(GridType& grid, const bool greenClosure = false )
{
  using namespace Dune;

  // skip empty grids
  if (grid.levelGridView(0).template begin<0>() == grid.levelGridView(0).template end<0>())
    return;

  // save start level
  const int startLevel = grid.maxLevel();
  // save start grid size
  const int startSize  = grid.size( 0 );

  // some things are different for bisection grids
  const bool bisectionGrid =
    Capabilities::isLevelwiseConforming<GridType> :: v  == false &&
    Capabilities::isLeafwiseConforming<GridType>  :: v  == true ;

  /*
     run the refinement check
   */
  checkAdaptRefinement(grid, greenClosure);

  /*
     run coarsening check
   */
  int counter = 0;
  const int counterEstimate = (grid.maxLevel() - startLevel) * 10;
  // now the same with coarsening
  while ( grid.maxLevel() > startLevel )
  {
    const int oldMaxLevel = grid.maxLevel();

    // mark all leaf entities for coarsening
    markForAdaptation( grid , -1 );

    bool markedCoarsen = grid.preAdapt();
    if( markedCoarsen == false )
    {
      DUNE_THROW(InvalidStateException,"grid.preAdapt() does not return correct information");
    }

    // check mightVanish
    typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
    ElementIterator eEndIt = grid.levelGridView(0).template end<0>();

    for (ElementIterator it = grid.levelGridView(0).template begin<0>();
         it != eEndIt; ++ it)
    {
      checkHierarchy( *it );
    }

    // only marked for coarsening ==> refined = false
    bool refined = grid.adapt() ;
    if( refined )
    {
      DUNE_THROW(InvalidStateException,"grid.adapt() returns wrong information");
    }

    /// check new max level
    if( grid.maxLevel () >= oldMaxLevel )
    {
      if( ! bisectionGrid )
        DUNE_THROW(InvalidStateException,"grid.maxLevel() wrong after coarsening " << oldMaxLevel << " " << grid.maxLevel() );
    }

    // Loop over all levels except the lowest one
    for (int level = 0 ; level <= grid.maxLevel(); ++level )
    {
      typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
      ElementIterator eEndIt = grid.levelGridView(level).template end<0>();

      for (ElementIterator it = grid.levelGridView(level).template begin<0>();
           it != eEndIt; ++ it)
      {
        if( it->isNew () )
        {
          DUNE_THROW(InvalidStateException,"After postAdapt() was called no entity is new, i.e., isNew() == false");
        }
      }
    }

    grid.postAdapt();

    // Loop over all levels except the lowest one
    for (int level = 0 ; level <= grid.maxLevel(); ++level )
    {
      typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
      ElementIterator eEndIt = grid.levelGridView(level).template end<0>();

      for (ElementIterator it = grid.levelGridView(level).template begin<0>();
           it != eEndIt; ++ it)
      {
        if( it->isNew () )
        {
          DUNE_THROW(InvalidStateException,"After postAdapt() was called no entity is new, i.e., isNew() == false");
        }
        if( it->mightVanish() )
        {
          DUNE_THROW(InvalidStateException,"After postAdapt() was called no entity might vanish, i.e., mightVanish() == false");
        }
      }
    }

    ++counter;
    if( counter > counterEstimate )
      DUNE_THROW(InvalidStateException,"Coarsening does not get back to startLevel");
  }

  const int newSize = grid.size( 0 );
  if( startSize != newSize )
  {
    dwarn << "After coarsening a different number of elements is obtained! old = "
          << startSize << "  new = " << newSize << std::endl;
  }

  if( startLevel != grid.maxLevel() )
  {
    dwarn << "After coarsening a different maxLevel is obtained! old = "
          << startLevel << "  new = " << grid.maxLevel() << std::endl;
  }
}

#endif
