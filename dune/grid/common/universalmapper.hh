// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_UNIVERSALMAPPER_HH
#define DUNE_GRID_COMMON_UNIVERSALMAPPER_HH

#include <iostream>
#include <map>

#include <dune/common/deprecated.hh>

#include "mapper.hh"

#warning "<dune/grid/common/universalmapper.hh> is deprecated in DUNE 2.6"

/**
 * @file
 * @brief  Mapper for any set of entities
 * @author Peter Bastian
 */

namespace Dune
{
  /**
   * @addtogroup Mapper
   *
   * @{
   */

  /** @brief Implements a mapper for an arbitrary subset of entities

      This implementation uses an ID set and a map, thus it has log complexity for each access.
          Template parameters are:

      Entities need to be registered in order to use them. If an entity is queried with map, the known index is returned or a new index is created. The method contains only return true, if the entites was queried via map already.

   * \tparam G   A Dune grid type.
   * \tparam IDS An Id set type for the given grid.
   * \tparam IndexType Number type used for the indices
   */
  template <typename G, typename IDS, typename IndexType=int>
  class UniversalMapper :
    public Mapper<G,UniversalMapper<G,IDS> >
  {
    typedef typename IDS::IdType IdType;
  public:

    /** \brief Number type used for indices */
    typedef IndexType Index;

    /** @brief Construct mapper from grid and one of its id sets

       \param grid A Dune grid object.
       \param idset An IndexSet object of the grid.

     */
    DUNE_DEPRECATED_MSG("UniversalMapper is deprecated in DUNE 2.6")
    UniversalMapper (const G& grid, const IDS& idset)
      : g(grid), ids(idset), index_()
    {
      n=0;     // zero data elements
    }

    /** @brief Map entity to array index.

       If an entity is queried with map, the known index is returned or a new index is created. A call to map can never fail.

            \param e Reference to codim cc entity, where cc is the template parameter of the function.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class EntityType>
    Index index (const EntityType& e) const
    {
      IdType id = ids.id(e);                                 // get id
      typename std::map<IdType,Index>::iterator it = index_.find(id);    // look up in map
      if (it!=index_.end()) return it->second;               // return index if found
      index_[id] = n++;                                      // put next index in map
      return n-1;                                            // and return it
    }

    /** @brief Map subentity of codim 0 entity to array index.

       If an entity is queried with map, the known index is returned or a new index is created. A call to map can never fail.

       \param e Reference to codim 0 entity.
       \param i Number of codim cc subentity of e, where cc is the template parameter of the function.
       \param cc codim of the subentity
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    Index subIndex (const typename G::Traits::template Codim<0>::Entity& e, int i, int cc) const
    {
      IdType id = ids.subId(e,i,cc);           // get id
      typename std::map<IdType,Index>::iterator it = index_.find(id);    // look up in map
      if (it!=index_.end()) return it->second;               // return index if found
      index_[id] = n++;                                      // put next index in map
      return n-1;                                            // and return it
    }

    /** @brief Return total number of entities in the entity set managed by the mapper.

       This number can be used to allocate a vector of data elements associated with the
       entities of the set. In the parallel case this number is per process (i.e. it
       may be different in different processes).

       \return Size of the entity set.
     */
    int size () const
    {
      return n;
    }

    /** @brief Returns true if the entity is contained in the index set

       The method contains only return true, if the entites was queried via map already.

       \param e Reference to entity
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<class EntityType>
    bool contains (const EntityType& e, Index& result) const
    {
      IdType id = ids.id(e);                                 // get id
      typename std::map<IdType,Index>::iterator it = index_.find(id);    // look up in map
      if (it!=index_.end())
      {
        result = it->second;
        return true;
      }
      else
        return false;
    }

    /** @brief Returns true if the entity is contained in the index set

       \param[in] e Reference to codim 0 entity
       \param[in] i subentity number
       \param[in] cc subentity codim
       \param[out] result integer reference where corresponding index is stored if true
       \return true if entity is in entity set of the mapper
     */
    bool contains (const typename G::Traits::template Codim<0>::Entity& e, int i, int cc, Index& result) const
    {
      IdType id = ids.subId(e,i,cc);           // get id
      typename std::map<IdType,Index>::iterator it = index_.find(id);    // look up in map
      if (it!=index_.end())
      {
        result = it->second;
        return true;
      }
      else
        return false;
    }

    /** @brief Recalculates map after mesh adaptation
     */
    void update ()
    {     // nothing to do here
    }

    // clear the mapper
    void clear ()
    {
      index_.clear();
      n = 0;
    }

  private:
    mutable int n;     // number of data elements required
    const G& g;
    const IDS& ids;
    mutable std::map<IdType,Index> index_;
  };




  /** @brief Universal mapper based on global ids

     Template parameters are:

   * \tparam G A Dune grid type.
   */
  template <typename G>
  class GlobalUniversalMapper : public UniversalMapper<G,typename G::Traits::GlobalIdSet>
  {
  public:
    /* @brief The constructor
       @param grid A reference to a grid.
     */
    DUNE_DEPRECATED_MSG("GlobalUniversalMapper is deprecated in DUNE 2.6")
    GlobalUniversalMapper (const G& grid)
      : UniversalMapper<G,typename G::Traits::GlobalIdSet>(grid,grid.globalIdSet())
    {}
  };

  /** @brief Universal mapper based on local ids

     Template parameters are:

     \tparam G A Dune grid type.
   */
  template <typename G>
  class LocalUniversalMapper : public UniversalMapper<G,typename G::Traits::LocalIdSet>
  {
  public:
    /**
     * \brief The constructor
     * \param grid A reference to a grid.
     */
    DUNE_DEPRECATED_MSG("LocalUniversalMapper is deprecated in DUNE 2.6")
    LocalUniversalMapper (const G& grid)
      : UniversalMapper<G,typename G::Traits::LocalIdSet>(grid,grid.localIdSet())
    {}
  };


  /** @} */
}
#endif
