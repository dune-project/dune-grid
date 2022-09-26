// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_PERSISTENTCONTAINER_HH
#define DUNE_ALBERTA_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>

#if HAVE_ALBERTA
#include <dune/grid/utility/persistentcontainervector.hh>

namespace Dune
{

  // PersistentContainer for AlbertaGrid
  // -----------------------------------

  template< int dim, int dimworld, class T >
  class PersistentContainer< AlbertaGrid< dim, dimworld >, T >
    : public PersistentContainerVector< AlbertaGrid< dim, dimworld >, typename AlbertaGrid< dim, dimworld >::HierarchicIndexSet, std::vector< T > >
  {
    typedef PersistentContainerVector< AlbertaGrid< dim, dimworld >, typename AlbertaGrid< dim, dimworld >::HierarchicIndexSet, std::vector< T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid.hierarchicIndexSet(), codim, value )
    {}
  };

} // end namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
