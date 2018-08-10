// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_SCSGMAPPER_HH
#define DUNE_GRID_COMMON_SCSGMAPPER_HH

#include <iostream>
#include "mapper.hh"

#include <dune/grid/common/grid.hh>

/**
 * @file
 * @brief  Mapper classes are used to attach data to a grid
 * @author Peter Bastian
 */

namespace Dune
{
  /**
   * @addtogroup Mapper
   *
   * @{
   */

  /** @brief Implementation class for a single codim and single geometry type mapper.
   *
   * In this implementation of a mapper the entity set used as domain for the map consists
   * of the entities of a given codimension c for all entities in the given index set. The index
   * set may only contain entities of a single geometry type, otherwise an exception is thrown. This
   * version is usually not used directly but is used to implement versions for leafwise and levelwise
   * entity sets.
   *
   * \tparam GV A Dune grid view type
   *
   * \tparam c A valid codimension
   */
  template <typename GV, int c>
  class SingleCodimSingleGeomTypeMapper :
    public Mapper<typename GV::Grid,SingleCodimSingleGeomTypeMapper<GV,c>, typename GV::IndexSet::IndexType >
  {
  public:

    /** \brief Number type used for indices */
    typedef typename GV::IndexSet::IndexType Index;

    /** \brief Number type used for the overall size (the return value of the 'size' method)
     *
     * The type used here is set to be the corresponding type used by the GridView's index set.
     */
    using size_type = decltype(std::declval<typename GV::IndexSet>().size(0));

    /** @brief Construct mapper from grid and one of its index sets.

       \param gridView A Dune GridView object.
     */
    SingleCodimSingleGeomTypeMapper (const GV& gridView)
     : is(gridView.indexSet())
    {
      // check that grid has only a single geometry type
      if (is.types(c).size() != 1)
        DUNE_THROW(GridError, "mapper treats only a single codim and a single geometry type");
    }

    /** @brief Map entity to array index.

            \param e Reference to codim cc entity, where cc is the template parameter of the function.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class EntityType>
    Index index (const EntityType& e) const
    {
      static_assert(EntityType::codimension == c, "Entity of wrong codim passed to SingleCodimSingleGeomTypeMapper");
      return is.index(e);
    }

    /** @brief Map subentity of codim 0 entity to array index.

       \param e Reference to codim 0 entity.
       \param i Number of the subentity of e, where cc is the template parameter of the function.
       \param codim Codimension of the subentity of e
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    Index subIndex (const typename GV::template Codim<0>::Entity& e,
               int i, unsigned int codim) const
    {
      if (codim != c)
        DUNE_THROW(GridError, "Id of wrong codim requested from SingleCodimSingleGeomTypeMapper");
      return is.subIndex(e,i,codim);
    }

    /** @brief Return total number of entities in the entity set managed by the mapper.

       This number can be used to allocate a vector of data elements associated with the
       entities of the set. In the parallel case this number is per process (i.e. it
       may be different in different processes).

       \return Size of the entity set.
     */
    size_type size () const
    {
      return is.size(c);
    }

    /** @brief Returns true if the entity is contained in the index set

       \param e Reference to entity
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<class EntityType>
    bool contains (const EntityType& e, Index& result) const
    {
      result = index(e);
      return true;
    }

    /** @brief Returns true if the entity is contained in the index set

       \param e Reference to codim 0 entity
       \param i subentity number
       \param cc subentity codim
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    bool contains (const typename GV::template Codim<0>::Entity& e, int i, int cc, Index& result) const
    {
      result = subIndex(e,i,cc);
      return true;
    }

    /** @brief Recalculates map after mesh adaptation
     */
    void update ()
    {     // nothing to do here
    }

  private:
    const typename GV::IndexSet& is;
  };

  /** @} */

  /**
   * @addtogroup Mapper
   *
   * @{
   */
  /** @brief Single codim and single geometry type mapper for leaf entities.


     This mapper uses all leaf entities of a certain codimension as its entity set. It is
     assumed (and checked) that the given grid contains only entities of a single geometry type.

     Template parameters are:

   * \tparam G A Dune grid type.
   * \tparam c A valid codimension.
   */
  template <typename G, int c>
  class LeafSingleCodimSingleGeomTypeMapper : public SingleCodimSingleGeomTypeMapper<typename G::LeafGridView,c> {
  public:
    /** \brief The constructor
     * \param grid A reference to a grid.
     */
    LeafSingleCodimSingleGeomTypeMapper (const G& grid)
      : SingleCodimSingleGeomTypeMapper<typename G::LeafGridView,c>(grid.leafGridView())
    {}
  };

  /** @brief Single codim and single geometry type mapper for entities of one level.


     This mapper uses all entities of a certain codimension on a given level as its entity set. It is
     assumed (and checked) that the given grid contains only entities of a single geometry type.

     Template parameters are:

   * \tparam G A Dune grid type.
   * \tparam c A valid codimension.
   */
  template <typename G, int c>
  class LevelSingleCodimSingleGeomTypeMapper : public SingleCodimSingleGeomTypeMapper<typename G::LevelGridView,c> {
  public:
    /* @brief The constructor
       @param grid A reference to a grid.
       @param level A valid level of the grid.
     */
    LevelSingleCodimSingleGeomTypeMapper (const G& grid, int level)
      : SingleCodimSingleGeomTypeMapper<typename G::LevelGridView,c>(grid.levelGridView(level))
    {}
  };

  /** @} */
}
#endif
