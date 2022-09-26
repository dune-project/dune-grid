// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_HH
#define DUNE_ONE_D_GRID_HH

#include <tuple>
#include <vector>
#include <list>

#include <dune/common/parallel/communication.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/type.hh>

/** \file
 * \brief The OneDGrid class
 */

#include "onedgrid/onedgridlist.hh"
#include "onedgrid/nulliteratorfactory.hh"
#include "onedgrid/onedgridentity.hh"
#include "onedgrid/onedgridentityseed.hh"
#include "onedgrid/onedgridintersections.hh"
#include "onedgrid/onedgridintersectioniterators.hh"
#include "onedgrid/onedgridleafiterator.hh"
#include "onedgrid/onedgridviews.hh"
#include "onedgrid/onedgridleveliterator.hh"
#include "onedgrid/onedgridhieriterator.hh"
#include "onedgrid/onedgridindexsets.hh"

namespace Dune {

  class OneDGrid;

  /** \brief The type used to for OneDGrid geometries

    If you ever want OneDGrid to use a different type for coordinates,
    you need to change the first argument of AxisAlignedCubeGeometry here.
  */
  template <int mydim, int coorddim, class GridImp>
  using OneDGridGeometry = AxisAlignedCubeGeometry<double, mydim, coorddim>;

  struct OneDGridFamily
  {
    typedef GridTraits<1,   // Grid dimension
                       1,   // Dimension of the physical space
                       Dune::OneDGrid,
        OneDGridGeometry,
        OneDGridEntity,
        OneDGridLevelIterator,
        OneDGridLeafIntersection,
        OneDGridLevelIntersection,
        OneDGridLeafIntersectionIterator,
        OneDGridLevelIntersectionIterator,
        OneDGridHierarchicIterator,
        OneDGridLeafIterator,
        OneDGridLevelIndexSet<const OneDGrid>,
        OneDGridLeafIndexSet<const OneDGrid>,
        OneDGridIdSet<const OneDGrid>,
        unsigned int,
        OneDGridIdSet<const OneDGrid>,
        unsigned int,
        Communication<No_Comm>,
        OneDGridLevelGridViewTraits,
        OneDGridLeafGridViewTraits,
        OneDGridEntitySeed>
    Traits;
  };

  //**********************************************************************
  //
  // --OneDGrid
  //
  //**********************************************************************

  /**
     \brief One-dimensional adaptive grid

     [<em> provides \ref Dune::Grid </em>]
     \ingroup GridImplementations
     \ingroup OneDGrid

     This implementation of the grid interface provides one-dimensional
     grids only. The OneDGrid can be nonuniform
     and provides local mesh refinement and coarsening.
   */
  class OneDGrid : public GridDefaultImplementation <1, 1,typename OneDGridGeometry<0,1,OneDGrid>::ctype, OneDGridFamily>
  {
    // Grid and world dimension are hardwired in this grid
    constexpr static int dim = 1;
    constexpr static int dimworld = 1;

    template <int , PartitionIteratorType, class >
    friend class OneDGridLevelIterator;

    friend class OneDGridHierarchicIterator<const OneDGrid>;

    template <int codim_, int dim_, class GridImp_>
    friend class OneDGridEntity;
    friend class OneDGridHierarchicIterator<OneDGrid>;
    friend class OneDGridLeafIntersection<const OneDGrid>;
    friend class OneDGridLevelIntersection<const OneDGrid>;
    friend class OneDGridLeafIntersectionIterator<const OneDGrid>;
    friend class OneDGridLevelIntersectionIterator<const OneDGrid>;

    friend class OneDGridLevelIndexSet<const OneDGrid>;
    friend class OneDGridLeafIndexSet<const OneDGrid>;
    friend class OneDGridIdSet<const OneDGrid>;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class OneDGridLeafIterator;

    friend class OneDGridLeafGridView<const OneDGrid>;
    friend class OneDGridLevelGridView<const OneDGrid>;

    template <class GridType_>
    friend class GridFactory;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    /** \brief Default constructor for the GridFactory */
    OneDGrid();

    // **********************************************************
    // The Interface Methods
    // **********************************************************

  public:

    /** \brief The type used to store coordinates
     */
    typedef typename OneDGridGeometry<0,1,OneDGrid>::ctype ctype;

    /** \brief GridFamily of OneDGrid */
    typedef OneDGridFamily GridFamily;

    //Provides the standard grid types
    typedef OneDGridFamily::Traits Traits;

    /** \brief Constructor with an explicit set of coordinates */
    OneDGrid(const std::vector<ctype>& coords);

    /** \brief Constructor for a uniform grid */
    OneDGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary);

    //! Destructor
    ~OneDGrid();

    /** \brief Return maximum level defined in this grid.

       Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {return entityImps_.size()-1;}

    /** \brief Create an Entity from an EntitySeed */
    template <typename Seed>
    static typename Traits::template Codim<Seed::codimension>::Entity
    entity(const Seed& seed)
    {
      const int codim = Seed::codimension;
      return typename Traits::template Codim<codim>::Entity(OneDGridEntity<codim,dim,const OneDGrid>(seed.impl().target()));
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      switch (codim)
      {
        case 0:
          return elements(level).size();
        case 1:
          return vertices(level).size();
        default:
          return 0;
      }
    }



    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return leafIndexSet().size(codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      // There is only one type for each codim
      return size(level,1-type.dim());
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }

    /** \brief Return the number of coarse grid boundary segments.

       For this grid implementation, the return value is always 2, because only connected domains
       are supported, and then the coarse grid boundary consists of two points.
     */
    size_t numBoundarySegments() const
    {
      return 2;
    }

    /** \brief Get the set of global ids */
    const Traits::GlobalIdSet& globalIdSet() const
    {
      return idSet_;
    }

    /** \brief Get the set of local ids */
    const Traits::LocalIdSet& localIdSet() const
    {
      return idSet_;
    }

    /** \brief Get an index set for the given level */
    const Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (! levelIndexSets_[level]) {
        levelIndexSets_[level] =
          new OneDGridLevelIndexSet<const OneDGrid>(*this, level);
        levelIndexSets_[level]->update();
      }

      return * levelIndexSets_[level];
    }

    /** \brief Get an index set for the leaf level */
    const Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** \brief Mark entity for refinement
     *
     * \param refCount if >0 mark for refinement, if <0 mark for coarsening
     * \param e Entity to the entity you want to mark
     *
     * \return True, if marking was successful
     */
    bool mark(int refCount, const Traits::Codim<0>::Entity& e );

    /** \brief return current adaptation marker of given entity

        \param e Entity to the entity you want to mark

        \return int current adaptation marker of entity pointer e
     */
    int getMark(const Traits::Codim<0>::Entity& e ) const;

    //! Does nothing except return true if some element has been marked for refinement
    bool preAdapt();

    //! Triggers the grid refinement process
    bool adapt();

    /** \brief Adaptation post-processing: Reset all adaptation state flags */
    void postAdapt();

    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    /** \brief The different forms of grid refinement supported by OneDGrid */
    enum RefinementType {
      /** \brief New level consists only of the refined elements */
      LOCAL,
      /** \brief New level consists of the refined elements and the unrefined ones, too */
      COPY
    };

    /** \brief Sets the type of grid refinement */
    void setRefinementType(RefinementType type) {
      refinementType_ = type;
    }

    /** \brief Does one uniform refinement step
     *
     * \param refCount I don't know what this is good for.  It doesn't
     *        actually do anything.
     */
    void globalRefine(int refCount);

    // dummy parallel functions

    const Communication &comm () const
    {
      return ccobj;
    }


  private:

    /** \brief Get vertex lists directly -- makes the code more readable */
    OneDGridList<OneDEntityImp<0> >& vertices(int level) {
      return std::get<0>(entityImps_[level]);
    }

    /** \brief Get vertex lists directly -- makes the code more readable */
    const OneDGridList<OneDEntityImp<0> >& vertices(int level) const {
      return std::get<0>(entityImps_[level]);
    }

    /** \brief Get element lists directly -- makes the code more readable */
    OneDGridList<OneDEntityImp<1> >& elements(int level) {
      return std::get<1>(entityImps_[level]);
    }

    /** \brief Get element lists directly -- makes the code more readable */
    const OneDGridList<OneDEntityImp<1> >& elements(int level) const {
      return std::get<1>(entityImps_[level]);
    }

    Communication ccobj;

    /** \brief Update all indices and ids */
    void setIndices();

    /** \brief Produce an entity id that has not been used in this grid before.
     */
    unsigned int getNextFreeId()
    {
      return freeIdCounter_++;
    }

    //! The type of grid refinement currently in use
    RefinementType refinementType_;

    OneDGridList<OneDEntityImp<0> >::iterator getLeftUpperVertex(const OneDEntityImp<1>* eIt);

    OneDGridList<OneDEntityImp<0> >::iterator getRightUpperVertex(const OneDEntityImp<1>* eIt);

    /** \brief Returns an iterator to the first element on the left of
        the input element which has sons.
     */
    OneDGridList<OneDEntityImp<1> >::iterator getLeftNeighborWithSon(OneDGridList<OneDEntityImp<1> >::iterator eIt);

    // The vertices and elements of the grid hierarchy
    std::vector<std::tuple<OneDGridList<OneDEntityImp<0> >,
            OneDGridList<OneDEntityImp<1> > > > entityImps_;

    // Our set of level indices
    mutable std::vector<OneDGridLevelIndexSet<const OneDGrid>* > levelIndexSets_;

    OneDGridLeafIndexSet<const OneDGrid> leafIndexSet_;

    OneDGridIdSet<const OneDGrid> idSet_;

    // Every entity gets a unique id, unless it is a copy of an entity on a coarser level.
    // This is the counter that we use to create the unique id.
    unsigned int freeIdCounter_;

    /** Since a OneDGrid is one-dimensional and connected, there can only be two possible numberings
        of the boundary segments.  Either the left one is '0' and the right one is '1' or the reverse.
        This flag stores which is the case. */
    bool reversedBoundarySegmentNumbering_;

  }; // end Class OneDGrid

  namespace Capabilities
  {
    /** \struct hasBackupRestoreFacilities
       \ingroup OneDGrid
     */

    /** \struct IsUnstructured
       \ingroup OneDGrid
     */

    /** \brief OneDGrid has only one geometry type for codim 0 entities
       \ingroup OneDGrid
     */
    template< >
    struct hasSingleGeometryType< OneDGrid >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(1).id();
    };


    /** \brief OneDGrid has entities for all codimension
       \ingroup OneDGrid
     */
    template<int cdim>
    struct hasEntity< OneDGrid, cdim >
    {
      static const bool v = true;
    };

    /**
     * \brief OneDGrid can iterate over all codimensions
     * \ingroup OneDGrid
     **/
    template<int codim>
    struct hasEntityIterator<OneDGrid, codim>
    {
      static const bool v = true;
    };

    /** \brief OneDGrid is levelwise conforming
       \ingroup OneDGrid
     */
    template<>
    struct isLevelwiseConforming< OneDGrid >
    {
      static const bool v = true;
    };

    /** \brief OneDGrid is leafwise conforming
       \ingroup OneDGrid
     */
    template<>
    struct isLeafwiseConforming< OneDGrid >
    {
      static const bool v = true;
    };

  }

} // namespace Dune

// Include the GridFactory specialization for OneDGrid, so everybody
// who includes the grid also gets the factory.  Since OneDGrid is
// not a template class, it needs to be a complete type before
// GridFactory<OneDGrid> can be defined.  This is why the #include-
// directive is at _the end_ of this file.
#include <dune/grid/onedgrid/onedgridfactory.hh>


#endif
