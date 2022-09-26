// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_UGGRID_HH
#define DUNE_UGGRID_HH

/** \file
 * \brief The UGGrid class
 */

#include <memory>

#include <dune/common/classname.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#if HAVE_DUNE_UGGRID || DOXYGEN

#ifdef ModelP
#include <dune/common/parallel/mpicommunication.hh>
#endif

/* [Before reading the following: the macros UG_DIM_2 and UG_DIM_3 where named
 *  _2 and _3, respectively, up until ug-3.12.0.]
 *
 * The following lines including the necessary UG headers are somewhat
   tricky.  Here's what's happening:
   UG can support two- and three-dimensional grids.  You choose by setting
   either UG_DIM_2 or UG_DIM_3 while compiling.  This changes all sorts of stuff, in
   particular data structures in the headers.
   UG was never supposed to provide 2d and 3d grids at the same time.
   However, when compiling it as c++, the dimension-dependent parts are
   wrapped up cleanly in the namespaces UG::D2 and UG::D3, respectively.  That
   way it is possible to link together the UG lib for 2d and the one for 3d.
   But we also need the headers twice!  Once with UG_DIM_2 set and once with UG_DIM_3!
   So here we go:*/

/* The following define tells the UG headers that we want access to a few
   special fields, for example the extra index fields in the element data structures.
   This define remains only for backwards compatibility with older version of UG.
   All dune-uggrid versions since 2016-08-05 do not need this #define or the #undef
   further below. */
#define FOR_DUNE

// Set UG's space-dimension flag to 2d
#define UG_DIM_2
// And include all necessary UG headers
#include "uggrid/ugincludes.hh"
#undef DUNE_UGINCLUDES_HH

// Wrap a few large UG macros by functions before they get undef'ed away.
// Here: The 2d-version of the macros
#define UG_DIM 2
#include "uggrid/ugwrapper.hh"
#undef UG_DIM

// UG defines a whole load of preprocessor macros.  ug_undefs.hh undefines
// them all, so we don't get name clashes.
#include "uggrid/ug_undefs.hh"
#undef UG_DIM_2

/* Now we're done with 2d, and we can do the whole thing over again for 3d */

/* All macros set by UG have been unset.  This includes the macros that ensure
   single inclusion of headers.  We can thus include them again.  However, we
   only want to include those headers again that contain dimension-dependent stuff.
   Therefore, we set a few single-inclusion defines manually before including
   ugincludes.hh again.
 */
#define UGTYPES_H
#define __HEAPS__
#define __UGENV__
#define __DEVICESH__
#ifdef ModelP
#define __PPIF__
#endif

#define UG_DIM_3
#include "uggrid/ugincludes.hh"
#undef DUNE_UGINCLUDES_HH

// Wrap a few large UG macros by functions before they get undef'ed away.
// This time it's the 3d-versions.
#define UG_DIM 3
#include "uggrid/ugwrapper.hh"
#undef UG_DIM

// undef all macros defined by UG
#include "uggrid/ug_undefs.hh"

#undef UG_DIM_3
#undef FOR_DUNE

// The components of the UGGrid interface
#include "uggrid/uggridgeometry.hh"
#include "uggrid/uggridlocalgeometry.hh"
#include "uggrid/uggridentity.hh"
#include "uggrid/uggridentityseed.hh"
#include "uggrid/uggridintersections.hh"
#include "uggrid/uggridintersectioniterators.hh"
#include "uggrid/uggridleveliterator.hh"
#include "uggrid/uggridleafiterator.hh"
#include "uggrid/uggridhieriterator.hh"
#include "uggrid/uggridindexsets.hh"
#include <dune/grid/uggrid/uggridviews.hh>
#ifdef ModelP
#include "uggrid/ugmessagebuffer.hh"
#include "uggrid/uglbgatherscatter.hh"
#endif

// Not needed here, but included for user convenience
#include "uggrid/uggridfactory.hh"

#ifdef ModelP
template <class DataHandle, int GridDim, int codim>
const Dune::UGGrid<GridDim>* Dune::UGMessageBuffer<DataHandle, GridDim, codim>::grid_;

template <class DataHandle, int GridDim, int codim>
DataHandle *Dune::UGMessageBuffer<DataHandle,GridDim,codim>::duneDataHandle_ = nullptr;

template <class DataHandle, int GridDim, int codim>
int Dune::UGMessageBuffer<DataHandle,GridDim,codim>::level = -1;
#endif // ModelP

namespace Dune {

#ifdef ModelP
    using UGCommunication = Communication<MPI_Comm>;
#else
    using UGCommunication = Communication<No_Comm>;
#endif

  template<int dim>
  struct UGGridFamily
  {
    typedef GridTraits<dim,dim,Dune::UGGrid<dim>,
        UGGridGeometry,
        UGGridEntity,
        UGGridLevelIterator,
        UGGridLeafIntersection,
        UGGridLevelIntersection,
        UGGridLeafIntersectionIterator,
        UGGridLevelIntersectionIterator,
        UGGridHierarchicIterator,
        UGGridLeafIterator,
        UGGridLevelIndexSet< const UGGrid<dim> >,
        UGGridLeafIndexSet< const UGGrid<dim> >,
        UGGridIdSet< const UGGrid<dim> >,
        typename UG_NS<dim>::UG_ID_TYPE,
        UGGridIdSet< const UGGrid<dim> >,
        typename UG_NS<dim>::UG_ID_TYPE,
        UGCommunication,
        UGGridLevelGridViewTraits,
        UGGridLeafGridViewTraits,
        UGGridEntitySeed,
        UGGridLocalGeometry>
    Traits;
  };


  //**********************************************************************
  //
  // --UGGrid
  //
  //**********************************************************************

  /**
     \brief Front-end for the grid manager of the finite element toolbox UG3

     \ingroup GridImplementations

     This is the implementation of the grid interface using the UG3 grid management
     system.  It is best described in this <a href="https://doi.org/10.1007/s007910050003">paper</a>.
     To our knowledge, the original code is not available anymore,
     but the relevant parts have been forked into the %Dune module
     dune-uggrid, available from <a href="https://www.dune-project.org/modules/dune-uggrid"></a>.

     UGGrid provides conforming grids
     in two and three space dimensions.  The grids can be mixed, i.e.
     2d grids can contain triangles and quadrilaterals and 3d grids can
     contain tetrahedra and hexahedra and also pyramids and prisms.
     The grid refinement rules are very flexible.  Local adaptive red/green
     refinement is the default, but a special method in the UGGrid class
     allows you to directly access a number of anisotropic refinements rules.
     Last but not least, the UG grid manager is completely parallelized,
     and you can use boundaries parametrized by either analytical expressions
     or high-resolution piecewise linear surfaces.

     In your %Dune application, you can instantiate objects of the
     type UGGrid<2> or UGGrid<3>.  You can have more than one, if
     you choose.  It is even possible to have 2d and 3d grids at the same
     time, even though the original UG system never intended to support
     this!

     See the documentation for the factory class GridFactory<UGGrid<dimworld> >
     to learn how to create UGGrid objects.
   */
  template <int dim>
  class UGGrid : public GridDefaultImplementation  <dim, dim, double, UGGridFamily<dim> >
  {
    typedef GridDefaultImplementation <dim, dim, double, UGGridFamily<dim> > Base;

    friend class UGGridGeometry<0,dim,const UGGrid<dim> >;
    friend class UGGridGeometry<dim,dim,const UGGrid<dim> >;
    friend class UGGridGeometry<1,2,const UGGrid<dim> >;
    friend class UGGridGeometry<2,3,const UGGrid<dim> >;

    friend class UGGridEntity <0,dim,const UGGrid<dim> >;
    friend class UGGridEntity <1,dim,const UGGrid<dim> >;
    friend class UGGridEntity <2,dim,const UGGrid<dim> >;
    friend class UGGridEntity <dim,dim,const UGGrid<dim> >;
    friend class UGGridHierarchicIterator<const UGGrid<dim> >;
    friend class UGGridLeafIntersection<const UGGrid<dim> >;
    friend class UGGridLevelIntersection<const UGGrid<dim> >;
    friend class UGGridLeafIntersectionIterator<const UGGrid<dim> >;
    friend class UGGridLevelIntersectionIterator<const UGGrid<dim> >;

    friend class UGGridLevelIndexSet<const UGGrid<dim> >;
    friend class UGGridLeafIndexSet<const UGGrid<dim> >;
    friend class UGGridIdSet<const UGGrid<dim> >;
    template <class GridImp_>
    friend class UGGridLeafGridView;
    template <class GridImp_>
    friend class UGGridLevelGridView;

    friend class GridFactory<UGGrid<dim> >;

#ifdef ModelP
    friend class UGLBGatherScatter;
#endif

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLeafIterator;
    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

    /** \brief UGGrid is only implemented for 2 and 3 dimension */
    static_assert(dim==2 || dim==3, "Use UGGrid only for 2d and 3d!");

    // The different instantiations are mutual friends so they can access
    // each others numOfUGGrids field
    friend class UGGrid<2>;
    friend class UGGrid<3>;

    //**********************************************************
    // The Interface Methods
    //**********************************************************
  public:
    //! type of the used GridFamily for this grid
    typedef UGGridFamily<dim>  GridFamily;

    // the Traits
    typedef typename UGGridFamily<dim>::Traits Traits;

    //! The type used to store coordinates
    typedef UG::DOUBLE ctype;

    //! The type used for process ranks
    typedef unsigned int Rank;

    /** \brief Default constructor
     */
    UGGrid(UGCommunication comm = {});

    //! Destructor
    ~UGGrid()  noexcept(false);

    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxlevel with 0 the coarsest level.
    int maxLevel() const;

    /** \brief Create an Entity from an EntitySeed */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::Entity
    entity(const Seed& seed) const
    {
      const int codim = Seed::codimension;
      return typename Traits::template Codim<codim>::Entity(UGGridEntity<codim,dim,const UGGrid<dim> >(seed.impl().target(),this));
    }

    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return leafIndexSet().size(codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return this->levelIndexSet(level).size(type);
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return this->leafIndexSet().size(type);
    }

    /** \brief Return the number of boundary segments */
    size_t numBoundarySegments() const {
      // The number is stored as a member of UGGrid upon grid creation.
      // The corresponding data structure is not exported by UG.  (It is in ug/dom/std/std_internal.h)
      return numBoundarySegments_;
    }

    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return idSet_;
    }

    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const
    {
      return idSet_;
    }

    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      return *levelIndexSets_[level];
    }

    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }

    /** @name Grid Refinement Methods */
    /*@{*/

    /** \brief Mark element for refinement
        \param refCount <ul>
        <li> 1: mark for red refinement </li>
        <li> -1: mark for coarsening </li>
        <li> 0: delete a possible refinement mark </li>
        </ul>
        \param e Element to be marked
       \return <ul>
        <li> true, if element was marked </li>
        <li> false, if nothing changed </li>
        </ul>
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e );

    /** \brief Mark method accepting a UG refinement rule

       \param e element to be marked for refinement
       \param rule One of the UG refinement rules
       \param side If rule==UG::%D2::%BLUE (one quadrilateral is split into two rectangles)
       you can choose the orientation of the cut by setting side==0 or side==1

       The available values for RefinementRule are:  (see the RefinementRule enum in ug/gm/gm.h)
       <h3>2D</h3>

       - NO_REFINEMENT
       - COPY
       - RED
       - BLUE
       - COARSE
       - BISECTION_1
       - BISECTION_2_Q
       - BISECTION_2_T1
       - BISECTION_2_T2
       - BISECTION_3

       <h3>3D</h3>

       - NO_REFINEMENT
       - COPY
       - RED
       - COARSE

       - TETRA_RED_HEX

       - PRISM_BISECT_1_2
       - PRISM_QUADSECT
       - PRISM_BISECT_HEX0
       - PRISM_BISECT_HEX1
       - PRISM_BISECT_HEX2
       - PRISM_ROTATE_LEFT
       - PRISM_ROTATE_RGHT
       - PRISM_QUADSECT_HEXPRI0
       - PRISM_RED_HEX
       - PRISM_BISECT_0_1
       - PRISM_BISECT_0_2
       - PRISM_BISECT_0_3

       - HEX_BISECT_0_1
       - HEX_BISECT_0_2
       - HEX_BISECT_0_3
       - HEX_TRISECT_0
       - HEX_TRISECT_5
       - HEX_QUADSECT_0
       - HEX_QUADSECT_1
       - HEX_QUADSECT_2
       - HEX_BISECT_HEXPRI0
       - HEX_BISECT_HEXPRI1

     */
    bool mark(const typename Traits::template Codim<0>::Entity & e,
              typename UG_NS<dim>::RefinementRule rule,
              int side=0);

    /** \brief Query whether element is marked for refinement */
    int getMark(const typename Traits::template Codim<0>::Entity& e) const;

    /** \brief returns true, if some elements might be coarsend during grid
       adaption, here always returns true */
    bool preAdapt();

    //! Triggers the grid refinement process
    bool adapt();

    /** \brief Clean up refinement markers */
    void postAdapt();
    /*@}*/

    /** \brief Distributes the grid and some data over the available nodes in a distributed machine

        \tparam DataHandle works like the data handle for the communicate
        methods.

        \return True, if grid has changed, false otherwise
     */
    template<class DataHandle>
    bool loadBalance (DataHandle& dataHandle)
    {
#ifdef ModelP
      // gather element data
      if (dataHandle.contains(dim, 0))
        UGLBGatherScatter::template gather<0>(this->leafGridView(), dataHandle);

      // gather node data
      if (dataHandle.contains(dim,dim))
        UGLBGatherScatter::template gather<dim>(this->leafGridView(), dataHandle);
#endif

      // the load balancing step now also attaches
      // the data to the entities and distributes it
      loadBalance();

#ifdef ModelP
      // scatter element data
      if (dataHandle.contains(dim, 0))
        UGLBGatherScatter::template scatter<0>(this->leafGridView(), dataHandle);

      // scatter node data
      if (dataHandle.contains(dim,dim))
        UGLBGatherScatter::template scatter<dim>(this->leafGridView(), dataHandle);
#endif

      return true;
    }

    /** \brief Distributes this grid over the available nodes in a distributed machine

       \bug The return value is always 'true'

       \param minlevel The coarsest grid level that gets distributed
     */
    bool loadBalance(int minlevel=0);

    /** \brief Distribute this grid over a distributed machine
     *
     * \param[in] targetProcessors For each leaf element the rank of the process the element shall be sent to
     * \param[in] fromLevel The lowest level that gets redistributed (set to 0 when in doubt)
     *
     * This method allows to (re-)distribute the grid controlled by an external grid repartitioning library.
     * You need to get that library to assign a target rank to each interior element in the leaf grid.  With this
     * information in a std::vector, call this method, and UG will do the actual repartitioning.
     * Each leaf element will be sent to the assigned target rank.  For all other elements we look at
     * where there children are being sent to.  The parent is then sent to where most of its children are
     * (Familienzusammenfuehrung).
     *
     * The size of the input array targetProcessors is expected to be equal to the number of elements in
     * the 'all'-partition, i.e., the number Interior elements plus the number of Ghost elements.
     * To get the array entry corresponding to an Interior element, a MultipleCodimMultipleGeomTypeMapper
     * with layout class MCMGElementLayout is used.  Array entries corresponding to Ghost elements are ignored.
     *
     * In some cases you may also want to leave the lowest levels on one process, to have them all together
     * for multigrid coarse grid corrections.  In that case, use the fromLevel parameter with a value other
     * than zero, to redistribute only elements above a certain level.
     *
     * The fromLevel argument is also needed to allow the compiler to distinguish this method from
     * the loadBalance method with a single template DataHandle argument.
     *
     * \note In theory you can assign a target rank to any element on any level, and UG will magically transfer
     * the element to that rank and make everything come out right.  This is not supported by the UGGrid interface,
     * because I didn't see a use case for it.  If you do need it please ask on the Dune mailing list.
     *
     * \return true
     */
    bool loadBalance(const std::vector<Rank>& targetProcessors, unsigned int fromLevel);

    /** \brief Distributes the grid over the processes of a parallel machine, and sends data along with it
     *
     * \param[in] targetProcessors For each leaf element the rank of the process the element shall be sent to
     * \param[in] fromLevel The lowest level that gets redistributed (set to 0 when in doubt)
     * \param[in,out] dataHandle A data handle object that does the gathering and scattering of data
     * \tparam DataHandle works like the data handle for the communicate methods.
     *
     * \return true
     */
    template<class DataHandle>
    bool loadBalance (const std::vector<Rank>& targetProcessors, unsigned int fromLevel, DataHandle& dataHandle)
    {
#ifdef ModelP
      // gather element data
      if (dataHandle.contains(dim, 0))
        UGLBGatherScatter::template gather<0>(this->leafGridView(), dataHandle);

      // gather node data
      if (dataHandle.contains(dim,dim))
        UGLBGatherScatter::template gather<dim>(this->leafGridView(), dataHandle);
#endif

      // the load balancing step now also attaches
      // the data to the entities and distributes it
      loadBalance(targetProcessors,fromLevel);

#ifdef ModelP
      // scatter element data
      if (dataHandle.contains(dim, 0))
        UGLBGatherScatter::template scatter<0>(this->leafGridView(), dataHandle);

      // scatter node data
      if (dataHandle.contains(dim,dim))
        UGLBGatherScatter::template scatter<dim>(this->leafGridView(), dataHandle);
#endif

      return true;
    }

    /** the communication */
    const UGCommunication& comm () const
    {
      return ccobj_;
    }

  protected:
#ifdef ModelP
    template <int codim, class GridView, class DataHandle>
    void communicateUG_(const GridView& gv, int level,
                        DataHandle &dataHandle,
                        InterfaceType iftype,
                        CommunicationDirection dir) const
    {
      typename UG_NS<dim>::DDD_IF_DIR ugIfDir;
      // Translate the communication direction from Dune-Speak to UG-Speak
      if (dir==ForwardCommunication)
        ugIfDir = UG_NS<dim>::IF_FORWARD();
      else
        ugIfDir = UG_NS<dim>::IF_BACKWARD();

      typedef UGMessageBuffer<DataHandle,dim,codim> UGMsgBuf;
      UGMsgBuf::duneDataHandle_ = &dataHandle;

      UGMsgBuf::level = level;

      std::vector<typename UG_NS<dim>::DDD_IF> ugIfs = findDDDInterfaces(iftype, codim);

      unsigned bufSize = UGMsgBuf::ugBufferSize(gv);
      if (!bufSize)
        return;     // we don't need to communicate if we don't have any data!
      UGMsgBuf::grid_ = this;
      for (unsigned i=0; i < ugIfs.size(); ++i)
        UG_NS<dim>::DDD_IFOneway(multigrid_->dddContext(),
                                 ugIfs[i],
                                 ugIfDir,
                                 bufSize,
                                 &UGMsgBuf::ugGather_,
                                 &UGMsgBuf::ugScatter_);
    }

    /** \brief Translate Dune communication interface to UG communication interface */
    std::vector<typename UG_NS<dim>::DDD_IF> findDDDInterfaces(InterfaceType iftype,
                                                               int codim) const;
#endif
  public:
    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    /** \brief Rudimentary substitute for a hierarchic iterator on faces
        \param e, elementSide Grid face specified by an element and one of its sides
        \param maxl The finest level that should be traversed by the iterator
        \param[out] childElements For each subface: element index, elementSide, and level
        \param[out] childElementSides Indices for transformation because Dune numbers the
                                      faces of several elements differently than UG
     */
    void getChildrenOfSubface(const typename Traits::template Codim<0>::Entity & e,
                              int elementSide,
                              int maxl,
                              std::vector<typename Traits::template Codim<0>::Entity>& childElements,
                              std::vector<unsigned char>& childElementSides) const;

    /** \brief The different forms of grid refinement that UG supports */
    enum RefinementType {
      /** \brief New level consists only of the refined elements and the closure*/
      LOCAL,
      /** \brief New level consists of the refined elements and the unrefined ones, too */
      COPY
    };

    /** \brief Decide whether to add a green closure to locally refined grid sections or not */
    enum ClosureType {
      /** \brief Standard red/green refinement */
      GREEN,
      /** \brief No closure, results in nonconforming meshes */
      NONE
    };

    /** \brief Sets the type of grid refinement */
    void setRefinementType(RefinementType type) {
      refinementType_ = type;
    }

    /** \brief Sets the type of grid refinement closure */
    void setClosureType(ClosureType type) {
      closureType_ = type;
    }

    /** \brief Sets a vertex to a new position

       Changing a vertex' position changes its position on all grid levels!*/
    void setPosition(const typename Traits::template Codim<dim>::Entity& e,
                     const FieldVector<double, dim>& pos);

    /** \brief Does uniform refinement
     *
     * \param n Number of uniform refinement steps
     */
    void globalRefine(int n);

    /** \brief Save entire grid hierarchy to disk

       Test implementation -- not working!
     */
    void saveState(const std::string& filename) const;

    /** \brief Read entire grid hierarchy from disk

       Test implementation -- not working!
     */
    void loadState(const std::string& filename);

  private:
    /** \brief UG multigrid, which contains the actual grid hierarchy structure */
    typename UG_NS<dim>::MultiGrid* multigrid_;

    /** \brief The communication object. */
    UGCommunication ccobj_;

    /** \brief Recomputes entity indices after the grid was changed
        \param setLevelZero If this is false, level indices of the level 0 are not touched
        \param nodePermutation Permutation array for the vertex level 0 indices.  If this is NULL,
        the identity is used.
     */
    void setIndices(bool setLevelZero,
                    std::vector<unsigned int>* nodePermutation);

    // Each UGGrid object has a unique name to identify it in the
    // UG environment structure
    std::string name_;

    // Our set of level indices
    std::vector<std::shared_ptr<UGGridLevelIndexSet<const UGGrid<dim> > > > levelIndexSets_;

    UGGridLeafIndexSet<const UGGrid<dim> > leafIndexSet_;

    // One id set implementation
    // Used for both the local and the global UGGrid id sets
    UGGridIdSet<const UGGrid<dim> > idSet_;

    //! The type of grid refinement currently in use
    RefinementType refinementType_;

    //! The type of grid refinement closure currently in use
    ClosureType closureType_;

    /** \brief Number of UGGrids currently in use.
     *
     * This counts the number of UGGrids currently instantiated.  All
     * constructors of UGGrid look at this variable.  If it is zero, they
     * initialize UG before proceeding.  Destructors use the same mechanism
     * to safely shut down UG after deleting the last UGGrid object.
     */
    static int numOfUGGrids;

    /** \brief Remember whether some element has been marked for refinement
        ever since the last call to adapt().

        This is here to implement the return value of adapt().
     */
    bool someElementHasBeenMarkedForRefinement_;

    /** \brief Remember whether some element has been marked for coarsening
        ever since the last call to adapt().

        This is here to implement the return value of preAdapt().
     */
    bool someElementHasBeenMarkedForCoarsening_;

    /** \brief The classes implementing the geometry of the boundary segments, if requested */
    std::vector<std::shared_ptr<BoundarySegment<dim> > > boundarySegments_;

    /** \brief Overall number of coarse grid boundary segments.

       This includes the number of linear segments.
       Hence numBoundarySegments_ >= boundarySegments_.size()   (greater than or equal)
     */
    unsigned int numBoundarySegments_;

  }; // end Class UGGrid

  namespace Capabilities
  {
    /** \struct hasEntity
       \ingroup UGGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup UGGrid
     */

    /** \struct IsUnstructured
       \ingroup UGGrid
     */

    /** \brief UGGrid has entities of all codimensions
       \ingroup UGGrid
     */
    template<int dim, int codim>
    struct hasEntity< UGGrid<dim>, codim>
    {
      static const bool v = true;
    };

    /**
     * \brief Set default for hasEntityIterator to false
     * UGGrid can currently only iterate over elements and vertices
     * \ingroup UGGrid
     **/
    template<int dim, int codim>
    struct hasEntityIterator<UGGrid<dim>, codim>
    {
      static const bool v = false;
    };

    /**
     * \brief UGGrid can iterate over codim=0 entities (elements)
     * \ingroup UGGrid
     **/
    template<int dim>
    struct hasEntityIterator<UGGrid<dim>, 0>
    {
      static const bool v = true;
    };

    /**
     * \brief UGGrid can iterate over codim=dim entities (vertices)
     * \ingroup UGGrid
     **/
    template<int dim>
    struct hasEntityIterator<UGGrid<dim>, dim>
    {
      static const bool v = true;
    };

    /** \brief UGGrid can communicate on entities of all (existing) codimensions
     *  \ingroup UGGrid
     */
    template<int dim, int codim>
    struct canCommunicate<UGGrid<dim>, codim>
    {
      static const bool v = (codim>=0 && codim<=dim);
    };

    /** \brief UGGrid is levelwise conforming
       \ingroup UGGrid
     */
    template<int dim>
    struct isLevelwiseConforming< UGGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief UGGrid may not be leafwise conforming
       \ingroup UGGrid
     */
    template<int dim>
    struct isLeafwiseConforming< UGGrid<dim> >
    {
      static const bool v = false;
    };

  }

} // namespace Dune

#endif   // HAVE_DUNE_UGGRID || DOXYGEN
#endif   // DUNE_UGGRID_HH
