// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINERMAP_HH
#define DUNE_PERSISTENTCONTAINERMAP_HH

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  // PersistentContainerMap
  // ----------------------

  /** \brief map-based implementation of the PersistentContainer */
  template< class G, class IdSet, class Map >
  class PersistentContainerMap
  {
    typedef PersistentContainerMap< G, IdSet, Map > This;

  protected:
    template< class reference, class iterator >
    class IteratorWrapper;

  public:
    typedef G Grid;

    typedef typename Map::mapped_type Value;
    typedef typename Map::size_type Size;

    typedef IteratorWrapper< const Value, typename Map::const_iterator > ConstIterator;
    typedef IteratorWrapper< Value, typename Map::iterator > Iterator;

    PersistentContainerMap ( const Grid &grid, int codim, const IdSet &idSet, const Value &value )
      : grid_( &grid ),
        codim_( codim ),
        idSet_( &idSet ),
        data_()
    {
      resize( value );
    }

    template< class Entity >
    const Value &operator[] ( const Entity &entity ) const
    {
      assert( Entity::codimension == codimension() );
      typename Map::const_iterator pos = data_.find( idSet().id( entity ) );
      assert( pos != data_.end() );
      return pos->second;
    }

    template< class Entity >
    Value &operator[] ( const Entity &entity )
    {
      assert( Entity::codimension == codimension() );
      typename Map::iterator pos = data_.find( idSet().id( entity ) );
      assert( pos != data_.end() );
      return pos->second;
    }

    template< class Entity >
    const Value &operator() ( const Entity &entity, int subEntity ) const
    {
      typename Map::const_iterator pos = data_.find( idSet().subId( entity, subEntity, codimension() ) );
      assert( pos != data_.end() );
      return pos->second;
    }

    template< class Entity >
    Value &operator() ( const Entity &entity, int subEntity )
    {
      typename Map::iterator pos = data_.find( idSet().subId( entity, subEntity, codimension() ) );
      assert( pos != data_.end() );
      return pos->second;
    }

    Size size () const { return data_.size(); }

    void resize ( const Value &value = Value() )
    {
      Hybrid::forEach( std::make_index_sequence< Grid::dimension+1 >{},
        [ & ]( auto i ){ if( int(i) == this->codimension() ) this->template resize< i >( value ); } );
    }

    void shrinkToFit () {}

    void fill ( const Value &value ) { std::fill( begin(), end(), value ); }

    void swap ( This &other )
    {
      std::swap( grid_, other.grid_ );
      std::swap( codim_, other.codim_ );
      std::swap( idSet_, other.idSet_ );
      std::swap( data_, other.data_ );
    }

    ConstIterator begin () const;
    Iterator begin ();

    ConstIterator end () const;
    Iterator end ();

    int codimension () const { return codim_; }

  protected:
    const Grid &grid () const { return *grid_; }

    template< int codim >
    void resize ( const Value &value );

    template< int codim >
    void migrateLevel ( int level, const Value &value, Map &data,
                        std::integral_constant< bool, true > );

    template< int codim >
    void migrateLevel ( int level, const Value &value, Map &data,
                        std::integral_constant< bool, false > );

    static void migrateEntry ( const typename IdSet::IdType &id, const Value &value,
                               Map &oldData, Map &newData );

    const IdSet &idSet () const { return *idSet_; }

    const Grid *grid_;
    int codim_;
    const IdSet *idSet_;
    Map data_;
  };



  // PersistentContainerMap::IteratorWrapper
  // ---------------------------------------

  template< class G, class IdSet, class Map >
  template< class value, class iterator >
  class PersistentContainerMap< G, IdSet, Map >::IteratorWrapper
    : public iterator
  {
    typedef IteratorWrapper< const value, typename Map::const_iterator > ConstWrapper;

  public:
    IteratorWrapper ( const iterator &it ) : it_( it ) {}

    operator ConstWrapper () const { return ConstWrapper( it_ ); }

    value &operator* () { return it_->second; }
    value *operator-> () { return &(it_->second); }

    bool operator== ( const IteratorWrapper &other ) const { return (it_ == other.it_); }
    bool operator!= ( const IteratorWrapper &other ) const { return (it_ != other.it_); }

    IteratorWrapper &operator++ () { ++it_; return *this; }

  private:
    iterator it_;
  };




  // Implementation of PersistentContainerMap
  // ----------------------------------------

  template< class G, class IdSet, class Map >
  inline typename PersistentContainerMap< G, IdSet, Map >::ConstIterator
  PersistentContainerMap< G, IdSet, Map >::begin () const
  {
    return ConstIterator( data_.begin() );
  }

  template< class G, class IdSet, class Map >
  inline typename PersistentContainerMap< G, IdSet, Map >::Iterator
  PersistentContainerMap< G, IdSet, Map >::begin ()
  {
    return Iterator( data_.begin() );
  }


  template< class G, class IdSet, class Map >
  inline typename PersistentContainerMap< G, IdSet, Map >::ConstIterator
  PersistentContainerMap< G, IdSet, Map >::end () const
  {
    return ConstIterator( data_.end() );
  }

  template< class G, class IdSet, class Map >
  inline typename PersistentContainerMap< G, IdSet, Map >::Iterator
  PersistentContainerMap< G, IdSet, Map >::end ()
  {
    return Iterator( data_.end() );
  }


  template< class G, class IdSet, class Map >
  template< int codim >
  inline void PersistentContainerMap< G, IdSet, Map >::resize ( const Value &value )
  {
    std::integral_constant< bool, Capabilities::hasEntityIterator< Grid, codim >::v > hasEntityIterator;
    assert( codim == codimension() );

    // create empty map and swap it with current map (no need to copy twice)
    Map data;
    std::swap( data, data_ );

    // copy all data from old map into new one (adding new entries, if necessary)
    const int maxLevel = grid().maxLevel();
    for ( int level = 0; level <= maxLevel; ++level )
      migrateLevel< codim >( level, value, data, hasEntityIterator );
  }


  template< class G, class IdSet, class Map >
  template< int codim >
  inline void PersistentContainerMap< G, IdSet, Map >
  ::migrateLevel ( int level, const Value &value, Map &data,
                   std::integral_constant< bool, true > )
  {
    typedef typename Grid::LevelGridView LevelView;
    typedef typename LevelView::template Codim< codim >::Iterator LevelIterator;

    const LevelView levelView = grid().levelGridView( level );
    const LevelIterator end = levelView.template end< codim >();
    for( LevelIterator it = levelView.template begin< codim >(); it != end; ++it )
      migrateEntry( idSet().id( *it ), value, data, data_ );
  }


  template< class G, class IdSet, class Map >
  template< int codim >
  inline void PersistentContainerMap< G, IdSet, Map >
  ::migrateLevel ( int level, const Value &value, Map &data,
                   std::integral_constant< bool, false > )
  {
    typedef typename Grid::LevelGridView LevelView;
    typedef typename LevelView::template Codim< 0 >::Iterator LevelIterator;

    const LevelView levelView = grid().levelGridView( level );
    const LevelIterator end = levelView.template end< 0 >();
    for( LevelIterator it = levelView.template begin< 0 >(); it != end; ++it )
    {
      const typename LevelIterator::Entity &entity = *it;
      const int subEntities = entity.subEntities( codim );
      for( int i = 0; i < subEntities; ++i )
        migrateEntry( idSet().subId( entity, i, codim ), value, data, data_ );
    }
  }


  template< class G, class IdSet, class Map >
  inline void PersistentContainerMap< G, IdSet, Map >
  ::migrateEntry ( const typename IdSet::IdType &id, const Value &value,
                   Map &oldData, Map &newData )
  {
    // insert entry for id
    const std::pair< typename Map::iterator, bool > inserted
      = newData.insert( std::make_pair( id, value ) );

    // if entry did not exist previously, copy data
    if( inserted.second )
    {
      const typename Map::iterator pos = oldData.find( id );
      if( pos != oldData.end() )
      {
        inserted.first->second = pos->second;
        oldData.erase( pos );
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERMAP_HH
