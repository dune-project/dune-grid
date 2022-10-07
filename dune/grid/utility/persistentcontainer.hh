// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <map>

#include <dune/grid/utility/persistentcontainermap.hh>

namespace Dune
{

  /** \brief A class for storing data during an adaptation cycle.
   *
   * \copydetails PersistentContainerInterface
   */
  template< class G, class T >
  class PersistentContainer
    : public PersistentContainerMap< G, typename G::LocalIdSet, std::map< typename G::LocalIdSet::IdType, T > >
  {
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::map< typename G::LocalIdSet::IdType, T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, grid.localIdSet(), value )
    {}
  };

  /** \brief refer PersistentContainer<const Grid> to the implementation of the mutable grid */
  template< class Grid, class T >
  class PersistentContainer< const Grid, T >
    : public PersistentContainer< Grid, T >
  {
    typedef PersistentContainer< Grid, T > Base;
  public:
    typedef typename Base::Value Value;

    PersistentContainer ( const typename Base::Grid &grid, int codim, const Value &value = Value() )
      : Base(grid, codim, value)
    {}
  };

} // namespace Dune


#if 0

// the following implementation can be used for a grid providing a hash for the id type

#include <unordered_map>

namespace Dune
{

  template< G, class T >
  class PersistentContainer
    : public PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > >
  {
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value )
      : Base( grid, codim, grid.localIdSet(), value )
    {}
  };

} // namespace Dune

#endif // #if 0

namespace std
{

  template< class G, class T >
  inline void swap ( Dune::PersistentContainer< G, T > &a, Dune::PersistentContainer< G, T > &b )
  {
    a.swap( b );
  }

} // namespace std

#endif // #ifndef DUNE_PERSISTENTCONTAINER_HH
