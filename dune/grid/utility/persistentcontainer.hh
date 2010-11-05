// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <vector>

namespace Dune {


  template <class Grid, int codim, class Data>
  class PersistentContainerVector
  {
  protected:
    typedef typename Grid :: HierarchicIndexSetType HierarchicIndexSetType;
    typedef Grid GridType;
    const GridType& grid_;
    const HierarchicIndexSetType& indexSet_;
    std::vector< Data > data_;

  public:
    typedef typename GridType :: template Codim< codim > :: Entity EntityType;
  public:
    PersistentContainer( const GridType & grid )
      : grid_( grid )
        , indexSet_( grid_.hierarchicIndexSet() )
        , data_()
    {
      // resize to current size
      adapt();
    }

    PersistentContainer( const PersistentContainer & other )
      : grid_( other.grid_ )
        , indexSet_( other.indexSet_ )
        , data_( other.data_ )
    {}

    Data& operator () (const EntityType& entity )
    {
      return data_[ indexSet_.index( entity ) ];
    }

    const Data& operator () (const EntityType& entity ) const
    {
      return data_[ indexSet_.index( entity ) ];
    }

    Data& operator () (const EntityType& entity, const int subEntity, const int codim )
    {
      return data_[ indexSet_.subIndex( entity, subEntity, codim ) ];
    }

    const Data& operator () (const EntityType& entity, const int subEntity, const int codim ) const
    {
      return data_[ indexSet_.subIndex( entity, subEntity, codim ) ];
    }

    void adapt() const
    {
      data_.resize( indexSet_.size( codim ) );
    }
  };

} // end namespace Dune
#endif // end DUNE_PERSISTENTCONTAINER_HH
