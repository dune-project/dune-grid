// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINERVECTOR_HH
#define DUNE_PERSISTENTCONTAINERVECTOR_HH

#include <algorithm>
#include <cassert>

namespace Dune
{

  // PersistentContainerVector
  // -------------------------

  /**
   * \brief vector-based implementation of the PersistentContainer
   *
   * Some grid implementations, like YaspGrid, can provide consecutive,
   * zero-based, persistent indices for the entire grid hierarchy.
   * This implementation of a PersistentContainer uses such an index set and
   * a std::vector-like container to store user data in an efficient and
   * persistent manner.
   *
   * \note The persistent index set is actually allowed to be non-consecutive,
   *       i.e., some indices might not be assigned to an entity.
   *       As the result of the size method on the index set is used to allocate
   *       storage for the vector, it must be larger than the largest used index.
   *
   * \note It is sufficient if the index set provides indices to the codimension
   *       the persistent container is created for.
   *       Neither the method types() nor the method contains() need to be
   *       implemented.
   *
   * \note The persistent container is currently restricted to index sets
   *       containing a single geometry type.
   *
   * \todo Actually, we use a mapper rather than an index set.
   *       This would automatically resolve two problems:
   *       - support multiple geometry types,
   *       - the requirement to store a reference to the index set
   *       .
   *
   * \tparam  G         type of grid
   * \tparam  IndexSet  type of persistent index set
   * \tparam  Vector    type of vector to store the data in
   */
  template< class G, class IndexSet, class Vector >
  class PersistentContainerVector
  {
    typedef PersistentContainerVector< G, IndexSet, Vector > This;

  public:
    typedef G Grid;

    typedef typename Vector::value_type Value;
    typedef typename Vector::size_type Size;
    typedef typename Vector::const_iterator ConstIterator;
    typedef typename Vector::iterator Iterator;

    typedef typename Vector::allocator_type Allocator;

    PersistentContainerVector ( const IndexSet &indexSet, int codim, const Value &value,
                                const Allocator &allocator = Allocator() )
      : codim_( codim ),
        indexSet_( &indexSet ),
        data_( indexSet.size( codim ), value, allocator )
    {}

    template< class Entity >
    const Value &operator[] ( const Entity &entity ) const
    {
      assert( Entity::codimension == codimension() );
      const Size index = indexSet().index( entity );
      assert( index < data_.size() );
      return data_[ index ];
    }

    template< class Entity >
    Value &operator[] ( const Entity &entity )
    {
      assert( Entity::codimension == codimension() );
      const Size index = indexSet().index( entity );
      assert( index < data_.size() );
      return data_[ index ];
    }

    template< class Entity >
    const Value &operator() ( const Entity &entity, int subEntity ) const
    {
      const Size index = indexSet().subIndex( entity, subEntity, codimension() );
      assert( index < data_.size() );
      return data_[ index ];
    }

    template< class Entity >
    Value &operator() ( const Entity &entity, int subEntity )
    {
      const Size index = indexSet().subIndex( entity, subEntity, codimension() );
      assert( index < data_.size() );
      return data_[ index ];
    }

    Size size () const { return data_.size(); }

    void resize ( const Value &value = Value() )
    {
      const Size indexSetSize = indexSet().size( codimension() );
      data_.resize( indexSetSize, value );
    }

    void shrinkToFit () {}

    void fill ( const Value &value ) { std::fill( begin(), end(), value ); }

    void swap ( This &other )
    {
      std::swap( codim_, other.codim_ );
      std::swap( indexSet_, other.indexSet_ );
      std::swap( data_, other.data_ );
    }

    ConstIterator begin () const { return data_.begin(); }
    Iterator begin () { return data_.begin(); }

    ConstIterator end () const { return data_.end(); }
    Iterator end () { return data_.end(); }

    int codimension () const { return codim_; }

  protected:
    const IndexSet &indexSet () const { return *indexSet_; }

    int codim_;
    const IndexSet *indexSet_;
    Vector data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERVECTOR_HH
