// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDINDEXSETS_HH
#define DUNE_ALU2DGRIDINDEXSETS_HH

//- System includes
#include <vector>

//- Dune includes
#include <dune/common/stdstreams.hh>
#include <dune/common/bigunsignedint.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/indexidset.hh>


//- Local includes
#include "alu2dinclude.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;

  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;



  // ALU2dGridHierarchicIndexSet
  // ---------------------------

  //! hierarchic index set of ALU2dGrid
  template <int dim, int dimworld, ALU2DSPACE ElementType eltype>
  class ALU2dGridHierarchicIndexSet :
    public IndexSet< ALU2dGrid< dim, dimworld, eltype >,
        ALU2dGridHierarchicIndexSet< dim, dimworld, eltype >, int >
  {
    typedef ALU2dGridHierarchicIndexSet< dim, dimworld, eltype > This;

    typedef ALU2dGrid< dim, dimworld, eltype > GridType;
    enum { numCodim = dim+1 }; // i.e. 3

    friend class ALU2dGrid< dim, dimworld, eltype >;

    ALU2dGridHierarchicIndexSet( const GridType &grid )
      : grid_( grid )
    {}

  public:
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

    //! return hierarchic index of given entity
    template< int codim >
    int index ( const typename GridType::Traits::template Codim< codim >::Entity &entity ) const
    {
      return GridType::getRealImplementation( entity ).getIndex();
    }

    //! return hierarchic index of given entity
    template< class Entity >
    int index ( const Entity &entity ) const
    {
      return GridType::getRealImplementation( entity ).getIndex();
    }

    //! return subIndex of given entity for codim sub entity
    int subIndex ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      return grid_.getRealImplementation( e ).subIndex( i, codim);
    }

    //! return size of indexset, i.e. maxindex+1
    //! for given type, if type is not exisiting within grid 0 is returned
    int size ( GeometryType type ) const
    {
      const int codim = dim-type.dim();
      assert( grid_.geomTypes(codim).size() == 1 );
      if( type != grid_.geomTypes(codim)[0] ) return 0;
      // return size of hierarchic index set
      return grid_.hierSetSize(codim);
    }

    //! return size of indexset, i.e. maxindex+1
    int size ( int codim ) const
    {
      // return size of hierarchic index set
      return grid_.hierSetSize(codim);
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_.geomTypes(codim);
    }

    //! return true because all entities are contained in this set
    template <class EntityType>
    bool contains (const EntityType &) const { return true; }

  private:
    // our Grid
    const GridType & grid_;
  };



  // ALU2dGridLocalIdSet
  // -------------------

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGridLocalIdSet
    : public IdSet< ALU2dGrid< dim, dimworld, eltype >, ALU2dGridLocalIdSet< dim, dimworld, eltype >, int >
  {
    typedef ALU2dGridLocalIdSet< dim, dimworld, eltype > This;
    typedef IdSet< ALU2dGrid< dim, dimworld, eltype >, This, int > Base;

    friend class ALU2dGrid< dim, dimworld, eltype >;

    typedef ALU2dGrid< dim, dimworld, eltype > Grid;

  public:
    //! export type of id
    typedef typename Base::IdType IdType;

  private:
    typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;

    // this means that only up to 536,870,912 entities are allowed
    static const IdType codimMultiplier = 1 << 29;

    // create local id set , only for the grid allowed
    ALU2dGridLocalIdSet ( const Grid &grid )
      : hset_( grid.hierarchicIndexSet() )
    {}

    // fake method to have the same method like GlobalIdSet
    void updateIdSet () {}

  public:

    //! return global id of given entity
    template< class Entity >
    IdType id ( const Entity &e ) const
    {
      return id< Entity::codimension >( e );
    }

    //! return global id of given entity
    template< int cd >
    IdType id ( const typename Grid::template Codim< cd >::Entity &e ) const
    {
      assert( hset_.size( cd ) < codimMultiplier );
      return cd*codimMultiplier + hset_.index( e );
    }

    //! return subId of given entity
    template< class Entity >
    IdType subId ( const Entity &e, int i, unsigned int codim ) const
    {
      return subId< Entity::codimension >( e, i, codim );
    }

    //! return subId of given entity
    template< int cd >
    IdType subId ( const typename Grid::template Codim< cd >::Entity &e, int i, unsigned int codim ) const
    {
      const int realCodim = cd+codim;
      assert( hset_.size( realCodim ) < codimMultiplier );
      return realCodim*codimMultiplier + hset_.subIndex( e, i, codim );
    }

  private:
    const HierarchicIndexSet &hset_;
  };

} // end namespace Dune

#endif // #ifndef DUNE_ALU2DGRIDINDEXSETS_HH
