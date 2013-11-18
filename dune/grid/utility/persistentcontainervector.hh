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

  /** \brief vector-based implementation of the PersistentContainer */
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


    // deprecated stuff, will be removed after Dune 2.3

    typedef Grid GridType DUNE_DEPRECATED_MSG("Use Grid instead.");
    typedef Value Data DUNE_DEPRECATED_MSG("Use Value instead.");

    void reserve () DUNE_DEPRECATED_MSG("Use resize() instead.")
    { return resize(); }

    void clear () DUNE_DEPRECATED_MSG("Use resize() instead.")
    {
      resize( Value() );
      shrinkToFit();
      fill( Value() );
    }

    void update () DUNE_DEPRECATED_MSG("Use resize() instead.")
    {
      resize( Value() );
      shrinkToFit();
    }

  protected:
    const IndexSet &indexSet () const { return *indexSet_; }

    int codim_;
    const IndexSet *indexSet_;
    Vector data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERVECTOR_HH
