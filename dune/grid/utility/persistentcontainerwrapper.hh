// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINERWRAPPER_HH
#define DUNE_PERSISTENTCONTAINERWRAPPER_HH

#include <dune/grid/utility/hostgridaccess.hh>
#include <dune/grid/utility/persistentcontainer.hh>

namespace Dune
{

  // PersistentContainerWrapper
  // --------------------------

  template< class G, class T >
  class PersistentContainerWrapper
  {
    typedef PersistentContainerWrapper< G, T > This;

    typedef Dune::HostGridAccess< G > HostGridAccess;

    typedef typename HostGridAccess::HostGrid HostGrid;
    typedef PersistentContainer< HostGrid, T > PersistentContainerHostGrid;

  public:
    typedef G Grid;

    typedef typename PersistentContainer< HostGrid, T >::Value Value;
    typedef typename PersistentContainer< HostGrid, T >::Size Size;

    typedef typename PersistentContainer< HostGrid, T >::Iterator Iterator;
    typedef typename PersistentContainer< HostGrid, T >::ConstIterator ConstIterator;

    PersistentContainerWrapper ( const Grid &grid, int codim, const Value &value = Value() )
      : hostContainer_( HostGridAccess::hostGrid( grid ), codim, value )
    {}

    template< class Entity >
    const Value &operator[] ( const Entity &entity ) const
    {
      return hostContainer_[ HostGridAccess::hostEntity( entity ) ];
    }

    template< class Entity >
    Value &operator[] ( const Entity &entity )
    {
      return hostContainer_[ HostGridAccess::hostEntity( entity ) ];
    }

    template< class Entity >
    const Value &operator() ( const Entity &entity, int subEntity ) const
    {
      return hostContainer_( HostGridAccess::hostEntity( entity ), subEntity );
    }

    template< class Entity >
    Value &operator() ( const Entity &entity, int subEntity )
    {
      return hostContainer_( HostGridAccess::hostEntity( entity ), subEntity );
    }

    Size size () const { return hostContainer_.size(); }

    void resize ( const Value &value = Value() ) { hostContainer_.resize( value ); }
    void shrinkToFit () { return hostContainer_.shrinkToFit(); }

    void fill ( const Value &value = Value() ) { hostContainer_.fill( value ); }

    void swap ( This &other ) { hostContainer_.swap( other.hostContainer_ ); }

    ConstIterator begin () const { return hostContainer_.begin(); }
    Iterator begin () { return hostContainer_.begin(); }

    ConstIterator end () const { return hostContainer_.end(); }
    Iterator end () { return hostContainer_.end(); }

    int codimension () const { return hostContainer_.codimension(); }

  protected:
    PersistentContainer< HostGrid, T > hostContainer_;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERWRAPPER_HH
