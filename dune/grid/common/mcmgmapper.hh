// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_COMMON_MCMGMAPPER_HH
#define DUNE_GRID_COMMON_MCMGMAPPER_HH

#include <functional>
#include <iostream>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
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

  //! Layout template for elements
  /**
   * This layout template is for use in the
   * MultipleCodimMultipleGeomTypeMapper.  It selects only elements (entities
   * with dim=dimgrid).
   *
   * \tparam dimgrid The dimension of the grid.
   *
   * \deprecated Use \ref mcmgElementLayout() instead.
   */
  template<int dimgrid> struct DUNE_DEPRECATED_MSG("The MCMG layout classes have been deprecated. Pass `mcmgElementLayout()` to the constructor instead")
  MCMGElementLayout {
    //! test whether entities of the given geometry type should be included in
    //! the map
    bool contains (Dune::GeometryType gt) const { return gt.dim()==dimgrid; }
  };

  //! Layout template for vertices
  /**
   * This layout template is for use in the
   * MultipleCodimMultipleGeomTypeMapper.  It selects only vertices (entities
   * with dim=0).
   *
   * \tparam dimgrid The dimension of the grid.
   *
   * \deprecated Use \ref mcmgVertexLayout() instead.
   */
  template<int dim> struct DUNE_DEPRECATED_MSG("The MCMG layout classes have been deprecated. Pass `mcmgVertexLayout()` to the constructor instead")
  MCMGVertexLayout {
    //! test whether entities of the given geometry type should be included in
    //! the map
    bool contains (Dune::GeometryType gt) const { return gt.dim()==0; }
  };

  namespace Impl {

    /*
     * Dummy layout to be used as the default for
     * `MultipleCodimMultipleGeomTypeMapper`.  It should never be used, but
     * we need a default.
     *
     * This class can be removed once the `LayoutClass` template parameter
     * of `MultipleCodimMultipleGeomTypeMapper` is removed.
     */
    template<int dimgrid>
    struct MCMGFailLayout {
      MCMGFailLayout()
        { DUNE_THROW(Exception, "The default layout class cannot be used"); }
      bool contains(GeometryType gt) const
        { DUNE_THROW(Exception, "The default layout class cannot be used"); }
    };

  } /* namespace Impl */

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
  using MCMGLayout = std::function<bool(GeometryType, int)>;

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
   * set may contain entities of several geometry types. This
   * version is usually not used directly but is used to implement versions for leafwise and levelwise
   * entity sets.
   *
   * The geometry types to be included in the mapper are selected using a
   * layout functional (\ref MCMGLayout) that is passed to the constructor.
   *
   * \tparam GV     A Dune GridView type.
   * \tparam LayoutClass (deprecated) A helper class template with a method contains(), that
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
   */
  template <typename GV, template<int> class LayoutClass = Impl::MCMGFailLayout>
  class MultipleCodimMultipleGeomTypeMapper :
    public Mapper<typename GV::Grid,MultipleCodimMultipleGeomTypeMapper<GV,LayoutClass>, typename GV::IndexSet::IndexType >
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
     *
     * \param gridView_ A Dune GridView object.
     * \param layout A layout object.
     *
     * \deprecated Use the constructor taking a \ref MCMGLayout instead.
     */
    MultipleCodimMultipleGeomTypeMapper(const GV& gridView_, const LayoutClass<GV::dimension> layout = {})
      DUNE_DEPRECATED_MSG("Use the constructor taking a `MCMGLayout` functional instead")
      : MultipleCodimMultipleGeomTypeMapper(gridView_, wrapLayoutClass(layout))
    {}

    /**
     * \brief construct mapper from grid and layout description
     *
     * The `layout` parameter is a functional describing entities of which
     * geometry types are included in the mapper.  For commonly used cases,
     * convenience functions are provided.  See the \ref MCMGLayout type
     * documentation for details.
     *
     * \param gridView grid view whose entities should be included in the mapper
     * \param layout   functional describing which entities are included in the mapper
     */
    MultipleCodimMultipleGeomTypeMapper(const GV& gridView, const MCMGLayout& layout)
      : gridView(gridView)
      , is(gridView.indexSet())
      , layout_(layout)
    {
      update();
    }

    /*!
     * \brief Map entity to array index.
     *
     * \tparam EntityType
     * \param e Reference to codim \a EntityType entity.
     * \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class EntityType>
    Index index (const EntityType& e) const
    {
      const GeometryType gt = e.type();
      assert(offset(gt) != invalidOffset);
      return is.index(e) + offset(gt);
    }

    /** @brief Map subentity of codim 0 entity to array index.

       \param e Reference to codim 0 entity.
       \param i Number of subentity of e
       \param codim Codimension of the subentity
       \return An index in the range 0 ... Max number of entities in set - 1.
     */
    Index subIndex (const typename GV::template Codim<0>::Entity& e, int i, unsigned int codim) const
    {
      const GeometryType eType = e.type();
      GeometryType gt = eType.isNone() ?
        GeometryType( GeometryType::none, GV::dimension - codim ) :
        ReferenceElements<double,GV::dimension>::general(eType).type(i,codim) ;
      //GeometryType gt=ReferenceElements<double,GV::dimension>::general(e.type()).type(i,codim);
      assert(offset(gt) != invalidOffset);
      return is.subIndex(e, i, codim) + offset(gt);
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

    /** @brief Returns true if the entity is contained in the index set

       \param e Reference to entity
       \param result integer reference where corresponding index is  stored if true
       \return true if entity is in entity set of the mapper
     */
    template<class EntityType>
    bool contains (const EntityType& e, Index& result) const
    {
      if(!is.contains(e) || offset(e.type()) == invalidOffset)
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
        GeometryType( GeometryType::none, GV::dimension - cc ) :
        ReferenceElements<double,GV::dimension>::general(eType).type(i,cc) ;
      //GeometryType gt=ReferenceElements<double,GV::dimension>::general(e.type()).type(i,cc);
      if (offset(gt) == invalidOffset)
        return false;
      result = is.subIndex(e, i, cc) + offset(gt);
      return true;
    }

    /** @brief Recalculates map after mesh adaptation
     */
    void update ()
    {
      n = 0;

      for (unsigned int codim = 0; codim <= GV::dimension; ++codim)
      {
        // walk over all geometry types in the codimension
        for (const GeometryType& gt : is.types(codim)) {
          Index offset;

          // if the geometry type is contained in the layout, increment offset
          if (layout()(gt, GV::Grid::dimension)) {
            offset = n;
            n += is.size(gt);
          }
          else {
            offset = invalidOffset;
          }

          offsets[GlobalGeometryTypeIndex::index(gt)] = offset;
        }
      }
    }

    const MCMGLayout &layout () const { return layout_; }

  private:
    Index offset(GeometryType gt) const
      { return offsets[GlobalGeometryTypeIndex::index(gt)]; }

    static const Index invalidOffset = std::numeric_limits<Index>::max();

    // number of data elements required
    unsigned int n;
    // GridView is needed to keep the IndexSet valid
    const GV gridView;
    const typename GV::IndexSet& is;
    // provide an array for the offsets
    std::array<Index, GlobalGeometryTypeIndex::size(GV::dimension)> offsets;
    const MCMGLayout layout_;     // get layout object

  protected:
    /**
     * \brief wrap legacy layout classes
     */
    static MCMGLayout wrapLayoutClass(const LayoutClass<GV::dimension>& layout)
    {
      /* `mutable` as the `contains()` method is not required to be const */
      return [layout = layout](GeometryType gt, int) mutable {
        return layout.contains(gt);
      };
    }
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  Leaf and level mapper
  //

  /** @brief Multiple codim and multiple geometry type mapper for leaf entities.

     This mapper uses all leaf entities of a certain codimension as its entity set.

     \tparam G      A %Dune grid type.
     \tparam LayoutClass (deprecated) A helper class template which determines which types of
                 entities are mapped by this mapper.  See
                 MultipleCodimMultipleGeomTypeMapper for how exactly this
                 template should look.
   */
  template <typename G, template<int> class LayoutClass = Impl::MCMGFailLayout>
  class LeafMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,LayoutClass>
  {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,
        LayoutClass> Base;
  public:
    /** @brief The constructor
     *
     * @param grid A reference to a grid.
     * @param layout A layout object
     *
     * \deprecated Use the constructor taking a \ref MCMGLayout instead.
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid, const LayoutClass<G::dimension> layout = {})
      DUNE_DEPRECATED_MSG("Use the constructor taking a `MCMGLayout` functional instead")
      : LeafMultipleCodimMultipleGeomTypeMapper(grid, Base::wrapLayoutClass(layout))
    {}

    /**
     * \brief constructor
     *
     * \param grid   reference to the grid
     * \param layout layout functional describing which geometry types to include in the map.
     */
    LeafMultipleCodimMultipleGeomTypeMapper (const G& grid, const MCMGLayout& layout)
      : Base(grid.leafGridView(), layout)
    {}
  };

  /** @brief Multiple codim and multiple geometry type mapper for entities of one level.


     This mapper uses all entities of a certain codimension on a given level as its entity set.

     \tparam G      A %Dune grid type.
     \tparam LayoutClass (deprecated) A helper class template which determines which types of
                 entities are mapped by this mapper.  See
                 MultipleCodimMultipleGeomTypeMapper for how exactly this
                 template should look.
   */
  template <typename G, template<int> class LayoutClass = Impl::MCMGFailLayout>
  class LevelMultipleCodimMultipleGeomTypeMapper
    : public MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,LayoutClass> {
    typedef MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,
        LayoutClass> Base;
  public:
    /** @brief The constructor
     *
     * @param grid A reference to a grid.
     * @param level A valid level of the grid.
     * @param layout A layout object
     *
     * \deprecated Use the constructor taking a \ref MCMGLayout instead.
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const LayoutClass<G::dimension> layout = {})
      DUNE_DEPRECATED_MSG("Use the constructor taking a `MCMGLayout` functional instead")
      : LevelMultipleCodimMultipleGeomTypeMapper(grid, level, Base::wrapLayoutClass(layout))
    {}

    /**
     * \brief constructor
     *
     * \param grid   reference to the grid
     * \param level  valid level of the grid
     * \param layout layout functional describing which geometry types to include in the map.
     */
    LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const MCMGLayout& layout)
      : Base(grid.levelGridView(level),layout)
    {}
  };

  /** @} */
}
#endif
