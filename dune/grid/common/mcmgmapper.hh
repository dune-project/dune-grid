// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MCMGMAPPER_HH
#define DUNE_MCMGMAPPER_HH

#include <iostream>
#include <map>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include "mapper.hh"

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

  //////////////////////////////////////////////////////////////////////
  //
  //  Common Layout templates
  //

  //! Layout template for elements
  /**
   * This layout template is for use in the
   * MultipleCodimMultipleGeomTypeMapper.  It selects only elements (entities
   * with dim=dimgrid).
   *
   * \tparam dimgrid The dimension of the grid.
   */
  template<int dimgrid> struct MCMGElementLayout {
    //! test whether entities of the given geometry type should be included in
    //! the map
    bool contains (Dune::GeometryType gt) { return gt.dim()==dimgrid; }
  };

  //! Layout template for vertices
  /**
   * This layout template is for use in the
   * MultipleCodimMultipleGeomTypeMapper.  It selects only vertices (entities
   * with dim=0).
   *
   * \tparam dimgrid The dimension of the grid.
   */
  template<int dim> struct MCMGVertexLayout {
    //! test whether entities of the given geometry type should be included in
    //! the map
    bool contains (Dune::GeometryType gt) { return gt.dim()==0; }
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  MultipleCodimMultipleGeomTypeMapper
  //

  /** @brief Implementation class for a multiple codim and multiple geometry type mapper.
   *
   * In this implementation of a mapper the entity set used as domain for the map consists
   * of the entities of a subset of codimensions in the given index set. The index
   * set may contain entities of several geometry types. This
   * version is usually not used directly but is used to implement versions for leafwise and levelwise
   * entity sets.
   *
   * \tparam GV     A Dune GridView type.
   * \tparam Layout A helper class template with a method contains(), that
   *                returns true for all geometry types that are in the domain
   *                of the map.  The class should be of the following shape
     \code
     template<int dimgrid>
     struct LayoutClass {
     bool contains (Dune::GeometryType gt) const {
        // Return true if gt is in the domain of the map
     }
     };
     \endcode
   *                The MultipleCodimMultipleGeomTypeMapper will always
   *                substitute the dimension of the grid for the template
   *                parameter dimgrid.
   *
   * If you don't want to use the default constructor of the LayoutClass you
   * can construct it yourself and hand it to the respective constructor (with
   * dimgrid=GV::%dimension).  In this case the layout class should be copy
   * constructible.
   *
   * There are two predefined Layout class templates for the common cases that
   * only elements or only vertices should be mapped: MCMGElementLayout and
   * MCMGVertexLayout.
   */
  template <typename GV, template<int> class Layout>
  class MultipleCodimMultipleGeomTypeMapper :
    public Mapper<typename GV::Grid,MultipleCodimMultipleGeomTypeMapper<GV,Layout> >
  {
  public:

    // the following lines need to be skipped for intel compilers, because they
    // lead to ambiguous calls to methods
#ifndef __INTEL_COMPILER
    //! import the base class implementation of map and contains (including the deprecated version)
    //! \todo remove in after next release
    using Mapper< typename GV::Grid, MultipleCodimMultipleGeomTypeMapper >::map;
    using Mapper< typename GV::Grid, MultipleCodimMultipleGeomTypeMapper >::contains;
#endif

    /** @brief Construct mapper from grid and one of its index sets.
     *
     * Use this constructor to provide a custom layout object e.g. not
     * using the default constructor.
     *
     * \param gridView A Dune GridView object.
     * \param layout A layout object.
     */
    MultipleCodimMultipleGeomTypeMapper (const GV& gridView, const Layout<GV::dimension> layout)
      : is(gridView.indexSet()), layout(layout)
    {
      update();
    }

    /** @brief Construct mapper from grid and one of its index sets.

       \param gridView A Dune GridView object.
     */
    MultipleCodimMultipleGeomTypeMapper (const GV& gridView)
      : is(gridView.indexSet())
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
       \param i Number of subentity of e
       \param codim Codimension of the subendity
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    int map (const typename GV::template Codim<0>::Entity& e, int i, unsigned int codim) const
    {
      GeometryType gt=ReferenceElements<double,GV::dimension>::general(e.type()).type(i,codim);
      std::map<GeometryType,int>::const_iterator it = offset.find(gt);
      assert(it!=offset.end());
      return is.subIndex(e,i,codim) + it->second;
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
       \param cc subentity codim
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    bool contains (const typename GV::template Codim<0>::Entity& e, int i, int cc, int& result) const
    {
      result = this->map(e,i,cc);
      return true;
    }


    /** @brief Recalculates map after mesh adaptation
     */
    void update ()
    {
      n=0;     // zero data elements
      for (int c=0; c<=GV::dimension; c++)
        offset.clear();         // clear all maps

      // Compute offsets for the different geometry types.
      // Note that mapper becomes invalid when the grid is modified.
      for (int c=0; c<=GV::dimension; c++)
        for (size_t i=0; i<is.geomTypes(c).size(); i++)
          if (layout.contains(is.geomTypes(c)[i]))
          {
            offset[is.geomTypes(c)[i]] = n;
            n += is.size(is.geomTypes(c)[i]);
          }
    }

  private:
    int n;     // number of data elements required
    const typename GV::IndexSet& is;
    std::map<GeometryType,int> offset;     // provide a map with all geometry types
    mutable Layout<GV::dimension> layout;     // get layout object
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  Leaf and level mapper
  //

  /** @brief Multiple codim and multiple geometry type mapper for leaf entities.

     This mapper uses all leaf entities of a certain codimension as its entity set.

     \tparam G      A %Dune grid type.
     \tparam Layout A helper class template which determines which types of
                 entities are mapped by this mapper.  See
                 MultipleCodimMultipleGeomTypeMapper for how exactly this
                 template should look.
   */
  template <typename G, template<int> class Layout>
  class LeafMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,Layout>
  {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,
        Layout> Base;
  public:
    /** @brief The constructor
         @param grid A reference to a grid.
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid)
      : Base(grid.leafView())
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
      : Base(grid.leafView(),layout)
    {}

  };

  /** @brief Multiple codim and multiple geometry type mapper for entities of one level.


     This mapper uses all entities of a certain codimension on a given level as its entity set.

     \tparam G      A %Dune grid type.
     \tparam Layout A helper class template which determines which types of
                 entities are mapped by this mapper.  See
                 MultipleCodimMultipleGeomTypeMapper for how exactly this
                 template should look.
   */
  template <typename G, template<int> class Layout>
  class LevelMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,Layout> {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,
        Layout> Base;
  public:
    /** @brief The constructor
         @param grid A reference to a grid.
         @param level A valid level of the grid.
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level)
      : Base(grid.levelView(level))
    {}

    /** @brief The constructor
     *
     * Use this constructor to provide a custom layout object e.g. not
     * using the default constructor.
     *
     * @param grid A reference to a grid.
     * @param level A valid level of the grid.
     * @param layout A layout object
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const Layout<G::dimension> layout)
      : Base(grid.levelView(level),layout)
    {}

  };

  /** @} */
}
#endif
