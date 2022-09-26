// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_COMMON_MCMGMAPPER_HH
#define DUNE_GRID_COMMON_MCMGMAPPER_HH

#include <functional>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/geometry/dimension.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

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

  /**
   * \brief layout function for `MultipleCodimMultipleGeomTypeMapper`
   *
   * The layout function indicates which entity types are contained in the
   * mapper.  It is called for each `GeometryType` and the grid dimension.
   * A `true` value indicates the geometry type is part of the map; `false`
   * that it is not.
   *
   * For commonly used layouts containing only entities of a fixed dimension
   * or codimension, convenience functions returning a `MCMGLayout` are provided.
   *
   * The following example is equivalent to \ref mcmgElementLayout():
     \code
     MCMGLayout layout = [](GeometryType gt, int griddim) {
       return gt.dim() == griddim;
     };
     \endcode
   *
   * \see MultipleCodimMultipleGeomTypeMapper
   * \see mcmgLayout(Codim<codim>)
   * \see mcmgLayout(Dim<dim>)
   * \see mcmgElementLayout()
   * \see mcmgVertexLayout()
   */
  using MCMGLayout = std::function<size_t(GeometryType, int)>;

  /**
   * \brief layout for entities of codimension `codim`
   *
   * \see MultipleCodimMultipleGeomTypeMapper
   */
  template<int codim>
  MCMGLayout mcmgLayout(Codim<codim>)
  {
    return [](GeometryType gt, int dimgrid) {
      return dimgrid - gt.dim() == codim;
    };
  }

  /**
   * \brief layout for entities of dimension `dim`
   *
   * \see MultipleCodimMultipleGeomTypeMapper
   */
  template<int dim>
  MCMGLayout mcmgLayout(Dim<dim>)
  {
    return [](GeometryType gt, int) {
      return gt.dim() == dim;
    };
  }

  /**
   * \brief layout for elements (codim-0 entities)
   *
   * \see MultipleCodimMultipleGeomTypeMapper
   */
  inline MCMGLayout mcmgElementLayout()
  {
    return mcmgLayout(Codim<0>());
  }

  /**
   * \brief layout for vertices (dim-0 entities)
   *
   * \see MultipleCodimMultipleGeomTypeMapper
   */
  inline MCMGLayout mcmgVertexLayout()
  {
    return mcmgLayout(Dim<0>());
  }

  //////////////////////////////////////////////////////////////////////
  //
  //  MultipleCodimMultipleGeomTypeMapper
  //

  /** @brief Implementation class for a multiple codim and multiple geometry type mapper.
   *
   * In this implementation of a mapper the entity set used as domain for the map consists
   * of the entities of a subset of codimensions in the given index set. The index
   * set may contain entities of several geometry types.
   *
   * The geometry types to be included in the mapper are selected using a
   * layout functional (\ref MCMGLayout) that is passed to the constructor.
   */
  template <typename GV>
  class MultipleCodimMultipleGeomTypeMapper :
    public Mapper<typename GV::Grid,MultipleCodimMultipleGeomTypeMapper<GV>, typename GV::IndexSet::IndexType >
  {
  public:

    /** \brief Underlying GridView */
    typedef GV GridView;

    /** \brief Number type used for indices */
    typedef typename GV::IndexSet::IndexType Index;

    /** \brief Number type used for the overall size (the return value of the 'size' method)
     *
     * The type used here is set to be the corresponding type used by the GridView's index set.
     */
    using size_type = decltype(std::declval<typename GV::IndexSet>().size(0));

    /**
     * \brief construct mapper from grid and layout description
     *
     * The `layout` parameter is a functional describing entities of which
     * geometry types are included in the mapper.  For commonly used cases,
     * convenience functions are provided.  See the \ref MCMGLayout type
     * documentation for details.
     *
     * \param gridView grid view whose entities should be included in the mapper
     * \param layout   functional describing how many dof to store on each entity (fixed per geometry type)
     */
    MultipleCodimMultipleGeomTypeMapper(const GV& gridView, const MCMGLayout& layout)
      : gridView_(gridView)
      , indexSet_(&gridView_.indexSet())
      , layout_(layout)
    {
      update(gridView);
    }

    /*!
     * \brief Map entity to starting index in array for dof block
     *
     * \tparam EntityType
     * \param e Reference to codim \a EntityType entity.
     * \return An index in the range 0 ... (Max number of entities in set)*blockSize - 1.
     */
    template<class EntityType>
    Index index (const EntityType& e) const
    {
      const GeometryType gt = e.type();
      assert(offset(gt) != invalidOffset);
      return indexSet_->index(e)*blockSize(gt) + offset(gt);
    }

    /** @brief Map subentity of codim 0 entity to starting index in array for dof block

       \param e Reference to codim 0 entity.
       \param i Number of subentity of e
       \param codim Codimension of the subentity
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    Index subIndex (const typename GV::template Codim<0>::Entity& e, int i, unsigned int codim) const
    {
      const GeometryType eType = e.type();
      GeometryType gt = eType.isNone() ?
        GeometryTypes::none( GV::dimension - codim ) :
        ReferenceElements<double,GV::dimension>::general(eType).type(i,codim) ;
      //GeometryType gt=ReferenceElements<double,GV::dimension>::general(e.type()).type(i,codim);
      assert(offset(gt) != invalidOffset);
      return indexSet_->subIndex(e, i, codim)*blockSize(gt) + offset(gt);
    }

    /** @brief Return total number of entities in the entity set managed by the mapper.

       This number can be used to allocate a vector of data elements associated with the
       entities of the set. In the parallel case this number is per process (i.e. it
       may be different in different processes).

       \return Size of the entity set.
     */
    size_type size () const
    {
      return n;
    }

    /** @brief return number of entries for a given geometry type **/
    size_type size(GeometryType gt) const
    {
      return blockSize(gt);
    }

    /** @brief return the geometry types with entries **/
    const std::vector< GeometryType >& types ( int codim ) const
    {
      return myTypes_[ codim ];
    }

    /** @brief Returns a pair with the starting point in the dof vector
     *         and the number of degrees of freedom if the entity is contained in the index set
     *         otherwise {0,0} is returned

       \param e Reference to entity
       \param result integer reference to the start of the block
       \return pair with first entry equal to index for that entity and the second entry
               the number of degrees of freedom (zero if entity is not in entity set of the mapper)
     */
    template<class EntityType>
    IntegralRange<Index> indices (const EntityType& e) const
    {
      if(!indexSet_->contains(e) || offset(e.type()) == invalidOffset)
        return {0,0};
      Index start = index(e);
      return {start, start+blockSize(e.type())};
    }

    /** @brief Returns a pair with the starting point in the dof vector
     *         and the number of degrees of freedom if the entity is contained in the index set
     *         otherwise {0,0} is returned

       \param e Reference to codim 0 entity
       \param i subentity number
       \param cc subentity codim
       \param result integer reference to the start of the block
       \return pair with first entry equal to index for that entity and the second entry
               the number of degrees of freedom (zero if sub entity is not in entity set of the mapper)
     */
    IntegralRange<Index> indices (const typename GV::template Codim<0>::Entity& e, int i, int cc) const
    {
      const GeometryType eType = e.type();
      const GeometryType gt = eType.isNone() ?
        GeometryTypes::none(GV::dimension - cc) :
        ReferenceElements<double,GV::dimension>::general(eType).type(i,cc) ;
      if (offset(gt) == invalidOffset)
        return {0,0};
      else
      {
        Index start = subIndex(e,i,cc);
        return {start, start+blockSize(gt)};
      }
    }

    /** @brief Returns true if the entity is contained in the index set

       \param e Reference to entity
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<class EntityType>
    bool contains (const EntityType& e, Index& result) const
    {
      if(!indexSet_->contains(e) || offset(e.type()) == invalidOffset)
      {
        result = 0;
        return false;
      }
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
      const GeometryType eType = e.type();
      const GeometryType gt = eType.isNone() ?
        GeometryTypes::none( GV::dimension - cc ) :
        ReferenceElements<double,GV::dimension>::general(eType).type(i,cc) ;
      if (offset(gt) == invalidOffset)
        return false;
      result = indexSet_->subIndex(e, i, cc)*blockSize(gt) + offset(gt);
      return true;
    }

    /** @brief Recalculates indices after grid adaptation
     *
     * After grid adaptation you need to call this to update
     * the stored gridview and recalculate the indices.
     */
    void update (const GV& gridView)
    {
      gridView_ = gridView;
      indexSet_ = &gridView_.indexSet();
      update_();
    }

    /** @brief Recalculates indices after grid adaptation
     *
     * After grid adaptation you need to call this to update
     * the stored gridview and recalculate the indices.
     */
    void update (GV&& gridView)
    {
      gridView_ = std::move(gridView);
      indexSet_ = &gridView_.indexSet();
      update_();
    }

    /** @brief Recalculates indices after grid adaptation
     */
    [[deprecated("Use update(gridView) instead! Will be removed after release 2.8.")]]
    void update ()
    {
      update_();
    }

    const MCMGLayout &layout () const { return layout_; }
    const GridView &gridView () const { return gridView_; }

  private:
    void update_()
    {
      n = 0;

      std::fill(offsets.begin(),offsets.end(),Index(0));
      std::fill(blocks.begin(),blocks.end(),Index(0));

      for (unsigned int codim = 0; codim <= GV::dimension; ++codim)
      {
        // walk over all geometry types in the codimension
        for (const GeometryType& gt : indexSet_->types(codim)) {
          Index offset;
          size_t block = layout()(gt, GV::Grid::dimension);

          // if the geometry type is contained in the layout, increment offset
          // and store geometry type
          if (block) {
            offset = n;
            n += indexSet_->size(gt) * block;
            myTypes_[codim].push_back(gt);
          }
          else {
            offset = invalidOffset;
          }

          offsets[GlobalGeometryTypeIndex::index(gt)] = offset;
          blocks[GlobalGeometryTypeIndex::index(gt)] = block;
        }
      }
    }

    Index offset(GeometryType gt) const
      { return offsets[GlobalGeometryTypeIndex::index(gt)]; }
    Index blockSize(GeometryType gt) const
      { return blocks[GlobalGeometryTypeIndex::index(gt)]; }

    static const Index invalidOffset = std::numeric_limits<Index>::max();

    // number of data elements required
    unsigned int n;
    // GridView is needed to keep the IndexSet valid
    GV gridView_;
    const typename GV::IndexSet* indexSet_;
    // provide an array for the offsets
    std::array<Index, GlobalGeometryTypeIndex::size(GV::dimension)> offsets;
    std::array<Index, GlobalGeometryTypeIndex::size(GV::dimension)> blocks;
    const MCMGLayout layout_;     // get layout object
    std::vector<GeometryType> myTypes_[GV::dimension+1];
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  Leaf and level mapper
  //

  /** @brief Multiple codim and multiple geometry type mapper for leaf entities.

     This mapper uses all leaf entities of a certain codimension as its entity set.
     \deprecated Use MultipleCodimMultipleGeomTypeMapper instead
     \tparam G      A %Dune grid type.
   */
  template <typename G>
  class [[deprecated("Use MultipleCodimMultipleGeomTypeMapper instead! Will be removed after release 2.8.")]]
  LeafMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView>
  {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView> Base;
  public:

    /**
     * \brief constructor
     *
     * \param grid   reference to the grid
     * \param layout layout functional describing which geometry types to include in the map.
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid, const MCMGLayout& layout)
      : Base(grid.leafGridView(), layout)
      , gridPtr_(&grid)
    {}

    /** @brief Recalculates indices after grid adaptation
     *
     * After grid adaptation you need to call this to update
     * the stored gridview and recalculate the indices.
     */
    void update ()
    {
      Base::update(gridPtr_->leafGridView());
    }

  private:
    const G* gridPtr_;
  };

  /** @brief Multiple codim and multiple geometry type mapper for entities of one level.


     This mapper uses all entities of a certain codimension on a given level as its entity set.
     \deprecated Use MultipleCodimMultipleGeomTypeMapper instead
     \tparam G      A %Dune grid type.
   */
  template <typename G>
  class [[deprecated("Use MultipleCodimMultipleGeomTypeMapper instead! Will be removed after release 2.8.")]]
  LevelMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView> {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView> Base;
  public:

    /**
     * \brief constructor
     *
     * \param grid   reference to the grid
     * \param level  valid level of the grid
     * \param layout layout functional describing which geometry types to include in the map.
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const MCMGLayout& layout)
      : Base(grid.levelGridView(level),layout)
      , gridPtr_(&grid)
      , level_(level)
    {}

    /** @brief Recalculates indices after grid adaptation
     *
     * After grid adaptation you need to call this to update
     * the stored gridview and recalculate the indices.
     */
    void update ()
    {
      Base::update(gridPtr_->levelGridView(level_));
    }

  private:
    const G* gridPtr_;
    int level_;
  };

  /** @} */
}
#endif
