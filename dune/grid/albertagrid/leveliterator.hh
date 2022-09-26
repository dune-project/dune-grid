// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_LEVELITERATOR_HH
#define DUNE_ALBERTA_LEVELITERATOR_HH

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/albertagrid/treeiterator.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // AlbertaGridLevelIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridLevelIterator
    : public AlbertaGridTreeIterator< codim, GridImp, false >
  {
    typedef AlbertaGridLevelIterator< codim, pitype, GridImp > This;
    typedef AlbertaGridTreeIterator< codim, GridImp, false > Base;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::MarkerVector MarkerVector;

    AlbertaGridLevelIterator ()
    {}

    //! Constructor making end iterator
    AlbertaGridLevelIterator ( const GridImp &grid, int level )
      : Base( grid, level )
    {}

    //! Constructor making begin iterator
    AlbertaGridLevelIterator ( const GridImp &grid,
                               const MarkerVector *vec,
                               int level )
      : Base( grid, vec, level )
    {}

    //! increment the iterator
    void increment ()
    {
      Base::increment();
    }
  };


  template< int codim, class GridImp >
  class AlbertaGridLevelIterator< codim, Ghost_Partition, GridImp >
    : public AlbertaGridTreeIterator< codim, GridImp, false >
  {
    typedef AlbertaGridLevelIterator< codim, Ghost_Partition, GridImp > This;
    typedef AlbertaGridTreeIterator< codim, GridImp, false > Base;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::MarkerVector MarkerVector;

    AlbertaGridLevelIterator ()
    {}

    //! Constructor making end iterator
    AlbertaGridLevelIterator ( const GridImp &grid, int level )
      : Base( grid, level )
    {}

    //! Constructor making begin iterator (which is the end iterator in this case)
    AlbertaGridLevelIterator ( const GridImp &grid,
                               const MarkerVector *,
                               int level )
      : Base( grid, level )
    {}

    //! increment the iterator
    void increment ()
    {
      Base::increment();
    }
  };
}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_LEVELITERATOR_HH
