// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDINDEXSET_HH
#define DUNE_GRID_YASPGRIDINDEXSET_HH

/** \file
 *
   \brief level-wise, non-persistent, consecutive indices for YaspGrid

 */

namespace Dune {

  template<class GridImp>
  class YaspIndexSet
    : public IndexSet< GridImp, YaspIndexSet< GridImp >, unsigned int >
  {
    typedef YaspIndexSet< GridImp > This;
    typedef IndexSet< GridImp, This, unsigned int > Base;

  public:
    typedef typename Base::IndexType IndexType;

    using Base::subIndex;

    //! constructor stores reference to a grid and level
    YaspIndexSet ( const GridImp &g, int l )
      : grid( g ),
        level( l )
    {
      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryType(GeometryType::cube,GridImp::dimension-codim));
    }

    //! get index of an entity
    template<int cc>
    IndexType index (const typename GridImp::Traits::template Codim<cc>::Entity& e) const
    {
      assert( cc == 0 || cc == GridImp::dimension );
      return grid.getRealImplementation(e).compressedIndex();
    }

    //! get index of subentity of an entity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      return grid.getRealImplementation(e).subCompressedIndex(i,codim);
    }

    //! get number of entities of given type and level (the level is known to the object)
    int size (GeometryType type) const
    {
      return grid.size( level, type );
    }

    //! return size of set for a given codim
    int size (int codim) const
    {
      return grid.size( level, codim );
    }

    //! return true if the given entity is contained in \f$E\f$.
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return e.level() == level;
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return mytypes[codim];
    }

  private:
    const GridImp& grid;
    int level;
    std::vector<GeometryType> mytypes[GridImp::dimension+1];
  };

}   // namespace Dune

#endif  // DUNE_GRID_YASPGRIDINDEXSET_HH
