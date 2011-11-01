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

  template< class G, class T, class Allocator >
  class PersistentContainerWrapper
  {
    typedef PersistentContainerWrapper< G, T, Allocator > This;

    typedef Dune::HostGridAccess< G > HostGridAccess;

    typedef typename HostGridAccess::HostGrid HostGrid;
    typedef PersistentContainer< HostGrid, T, Allocator > PersistentContainerHostGrid;

  public:
    typedef G Grid;
    typedef T Data;

    typedef typename PersistentContainerHostGrid::Iterator Iterator;
    typedef typename PersistentContainerHostGrid::ConstIterator ConstIterator;

    PersistentContainerWrapper ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
      : hostContainer_( HostGridAccess::hostGrid( grid ), codim, allocator )
    {}

    template< class Entity >
    Data &operator[] ( const Entity &entity )
    {
      return hostContainer_[ HostGridAccess::hostEntity( entity ) ];
    }

    template< class Entity >
    const Data &operator[] ( const Entity &entity ) const
    {
      return hostContainer_[ HostGridAccess::hostEntity( entity ) ];
    }

    template< class Entity >
    Data &operator() ( const Entity &entity, const int subEntity )
    {
      return hostContainer_( HostGridAccess::hostEntity( entity ), subEntity );
    }

    template< class Entity >
    const Data &operator() ( const Entity &entity, const int subEntity ) const
    {
      return hostContainer_( HostGridAccess::hostEntity( entity ), subEntity );
    }

    Iterator begin () { return hostContainer_.begin(); }
    ConstIterator begin () const { return hostContainer_.begin(); }

    Iterator end () { return hostContainer_.end(); }
    ConstIterator end () const { return hostContainer_.end(); }

    size_t size () const { return hostContainer_.size(); }

    void clear () { hostContainer_.clear(); }
    void reserve () { hostContainer_.reserve(); }
    void update () { hostContainer_.update(); }

  private:
    PersistentContainerHostGrid hostContainer_ ;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERWRAPPER_HH
