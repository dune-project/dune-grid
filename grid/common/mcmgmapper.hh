// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#ifndef DUNE_MCMGMAPPER_HH
#define DUNE_MCMGMAPPER_HH

#include <iostream>
#include <map>
#include "mapper.hh"
#include "referenceelements.hh"

/**
 * @file
 * @brief  Mapper for multiple codim and multiple geometry types
 * @author Peter Bastian
 */

namespace Dune
{
  /**
   * @addtogroup Mapper
   *
   * @{
   */

  /** @brief Implementation class for a multiple codim and multiple geometry type mapper.
   *
   * In this implementation of a mapper the entity set used as domain for the map consists
   * of the entities of a subset of codimensions in the given index set. The index
   * set may contain entities of several geometry types. This
   * version is usually not used directly but is used to implement versions for leafwise and levelwise
   * entity sets.
   *
   * Template parameters are:
   *
   * \par G
   *    A Dune grid type.
   * \par IS
   *    LeafIndexSet or LevelIndexSet type of the given grid
   * \par Layout
   *  A helper class with a method contains(), that returns true for all geometry
   *  types that are in the domain of the map.  The class should be of the following
   *  shape
     \code
     template<int dim>
     struct LayoutClass {
        bool contains (Dune::GeometryType gt) {
            // Return true if gt is in the domain of the map
        }
     };
     \endcode
   *
   * If you don't want to use the default constructor of the LayoutClass you can construct it yourself
   * and hand it to the respective constructor.
   */
  template <typename G, typename IS, template<int> class Layout>
  class MultipleCodimMultipleGeomTypeMapper : Mapper<G,MultipleCodimMultipleGeomTypeMapper<G,IS,Layout> > {
  public:

    /** @brief Construct mapper from grid and one of its index sets.
     *
     * Use this constructor to provide a custom layout object e.g. not
     * using the default constructor.
     *
     * \param grid A Dune grid object.
     * \param indexset IndexSet object returned by grid.
     * \param layout A layout object.
     */
    MultipleCodimMultipleGeomTypeMapper (const G& grid, const IS& indexset, const Layout<G::dimension> layout)
      : g(grid), is(indexset), layout(layout)
    {
      update();
    }

    /** @brief Construct mapper from grid and one of its index sets.

       \param grid A Dune grid object.
       \param indexset IndexSet object returned by grid.

     */
    MultipleCodimMultipleGeomTypeMapper (const G& grid, const IS& indexset)
      : g(grid), is(indexset)
    {
      update();
    }

    /** @brief Map entity to array index.

            \param e Reference to codim cc entity, where cc is the template parameter of the function.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class EntityType>
    int map (const EntityType& e) const
    {
      return is.index(e) + offset.find(e.type())->second;
    }


    /** @brief Map subentity of codim 0 entity to array index.

       \param e Reference to codim 0 entity.
       \param i Number of codim cc subentity of e, where cc is the template parameter of the function.
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<int cc>
    int map (const typename G::Traits::template Codim<0>::Entity& e, int i) const
    {
      GeometryType gt=ReferenceElements<double,G::dimension>::general(e.type()).type(i,cc);
      //	  std::cout << "map: cc=" << cc << " gt=" << gt << " offset=" << offset.find(gt)->second << std::endl;
      return is.template subIndex<cc>(e,i) + offset.find(gt)->second;
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

       \param e Reference to entity
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<class EntityType>
    bool contains (const EntityType& e, int& result) const
    {
      if(!is.contains(e) || !layout.contains(e.type()))
      {
        result = 0;
        return false;
      }
      result = map(e);
      return true;
    }

    /** @brief Returns true if the entity is contained in the index set

       \param e Reference to codim 0 entity
       \param i subentity number
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<int cc>     // this is now the subentity's codim
    bool contains (const typename G::Traits::template Codim<0>::Entity& e, int i, int& result) const
    {
      result = this->template map<cc>(e,i);
      return true;
    }


    /** @brief Recalculates map after mesh adaptation
     */
    void update ()
    {
      n=0;     // zero data elements
      for (int c=0; c<=G::dimension; c++)
        offset.clear();         // clear all maps

      // Compute offsets for the different geometry types.
      // Note that mapper becomes invalid when the grid is modified.
      for (int c=0; c<=G::dimension; c++)
        for (size_t i=0; i<is.geomTypes(c).size(); i++)
          if (layout.contains(is.geomTypes(c)[i]))
          {
            //			  std::cout << "offset " << c << " " << is.geomTypes(c)[i] << " is " << n << std::endl;
            offset[is.geomTypes(c)[i]] = n;
            n += is.size(is.geomTypes(c)[i]);
          }
    }

  private:
    int n;     // number of data elements required
    const G& g;
    const IS& is;
    std::map<GeometryType,int> offset;     // provide a map with all geometry types
    mutable Layout<G::dimension> layout;     // get layout object
  };

  /** @brief Multiple codim and multiple geometry type mapper for leaf entities.

     This mapper uses all leaf entities of a certain codimension as its entity set.

     Template parameters are:

     \par G
     A %Dune grid type.
     \par Layout
     A helper class with a method contains(), that returns true for all geometry
     types that are in the domain of the map.  The class should be of the following
     shape
     \code
     template<int dim>
     struct LayoutClass {
        bool contains (Dune::GeometryType gt) {
            // Return true if gt is in the domain of the map
        }
     };
     \endcode
   */
  template <typename G, template<int> class Layout>
  class LeafMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LeafIndexSet,Layout>
  {
  public:
    /** @brief The constructor
         @param grid A reference to a grid.
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid)
      : MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LeafIndexSet,Layout>(grid,grid.leafIndexSet())
    {}

    /** @brief The constructor
     *
     * Use this constructor to provide a custom layout object e.g. not
     * using the default constructor.
     *
     * @param grid A reference to a grid.
     * @param layout A layout object
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid, const Layout<G::dimension> layout)
      : MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LeafIndexSet,Layout>(grid,grid.leafIndexSet(),layout)
    {}

  };

  /** @brief Multiple codim and multiple geometry type mapper for entities of one level.


     This mapper uses all entities of a certain codimension on a given level as its entity set.

     Template parameters are:

     \par G
     A %Dune grid type.
     \par Layout
     A helper class with a method contains(), that returns true for all geometry
     types that are in the domain of the map.  The class should be of the following
     shape
     \code
     template<int dim>
     struct LayoutClass {
        bool contains (Dune::GeometryType gt) {
            // Return true if gt is in the domain of the map
        }
     };
     \endcode
   */
  template <typename G, template<int> class Layout>
  class LevelMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LevelIndexSet,Layout> {
  public:
    /** @brief The constructor
         @param grid A reference to a grid.
         @param level A valid level of the grid.
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level)
      : MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LevelIndexSet,Layout>(grid,grid.levelIndexSet(level))
    {}

    /** @brief The constructor
     *
     * Use this constructor to provide a custom layout object e.g. not
     * using the default constructor.
     *
     * @param grid A reference to a grid.
     * @param layout A layout object
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const Layout<G::dimension> layout)
      : MultipleCodimMultipleGeomTypeMapper<G,typename G::Traits::LevekIndexSet,Layout>(grid,grid.levelIndexSet(level),layout)
    {}

  };

  /** @} */
}
#endif
