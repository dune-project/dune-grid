// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDINDEXSET_HH
#define DUNE_GRID_YASPGRIDINDEXSET_HH

/** \file
 *
   \brief level-wise, non-persistent, consecutive indices for YaspGrid

 */

namespace Dune {

  /** \brief Implementation of Level- and LeafIndexSets for YaspGrid
   *
   * \tparam GridImp The YaspGrid class we are an index set for
   * \tparam isLeafIndexSet false: class functions as level index set,
   *         true: class functions as leaf index set
   */
  template<class GridImp, bool isLeafIndexSet>
  class YaspIndexSet
    : public IndexSet< GridImp, YaspIndexSet< GridImp, isLeafIndexSet >, unsigned int >
  {
    typedef YaspIndexSet< GridImp, isLeafIndexSet > This;
    typedef IndexSet< GridImp, This, unsigned int > Base;

  public:
    typedef typename Base::IndexType IndexType;

    using Base::subIndex;

    /** \brief Level grid view constructor stores reference to a grid and level */
    YaspIndexSet ( const GridImp &g, int l )
      : grid( g ),
        level( l )
    {
      assert(not isLeafIndexSet);

      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryTypes::cube(GridImp::dimension-codim));
    }

    /** \brief Level grid view constructor stores reference to a grid and level */
    YaspIndexSet ( const GridImp &g )
      : grid( g )
    {
      assert(isLeafIndexSet);

      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryTypes::cube(GridImp::dimension-codim));
    }

    //! get index of an entity
    template<int cc>
    IndexType index (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const
    {
      return e.impl().compressedIndex();
    }

    //! get index of subentity of an entity
    template< int cc >
    IndexType subIndex ( const typename std::remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      return e.impl().subCompressedIndex(i, codim);
    }

    //! get number of entities of given type and level (the level is known to the object)
    std::size_t size (GeometryType type) const
    {
      return (isLeafIndexSet)
        ? grid.size( type )
        : grid.size( level, type );
    }

    //! return size of set for a given codim
    std::size_t size (int codim) const
    {
      return (isLeafIndexSet)
        ? grid.size( codim )
        : grid.size( level, codim );
    }

    //! return true if the given entity is contained in \f$E\f$.
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return (isLeafIndexSet)
        ? e.level() == grid.maxLevel()
        : e.level() == level;
    }

    std::vector< GeometryType > types ( int codim ) const { return mytypes[ codim ]; }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return mytypes[codim];
    }

  private:
    const GridImp& grid;
    int level;
    std::vector<GeometryType> mytypes[std::remove_const<GridImp>::type::dimension+1];
  };

}   // namespace Dune

#endif  // DUNE_GRID_YASPGRIDINDEXSET_HH
