// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_UGGRID_HH
#define DUNE_UGGRID_HH

/** \file
 * \brief The UGGrid class
 */

#include <dune/common/classname.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#if HAVE_UG

#ifdef ModelP
#include <dune/common/mpicollectivecommunication.hh>
#endif

/* The following lines including the necessary UG headers are somewhat
   tricky.  Here's what's happening:
   UG can support two- and three-dimensional grids.  You choose be setting
   either _2 oder _3 while compiling.  This changes all sorts of stuff, in
   particular data structures in the headers.
   UG was never supposed to provide 2d and 3d grids at the same time.
   However, when compiling it as c++, the dimension-dependent parts are
   wrapped up cleanly in the namespaces UG::D2 and UG::D3, respectively.  That
   way it is possible to link together the UG lib for 2d and the one for 3d.
   But we also need the headers twice!  Once with _2 set and once with _3!
   So here we go:*/

/* The following define tells the UG headers that we want access to a few
   special fields, for example the extra index fields in the element data structures. */
#define FOR_DUNE

// Set UG's space-dimension flag to 2d
#define _2
// And include all necessary UG headers
#include "uggrid/ugincludes.hh"

// Wrap a few large UG macros by functions before they get undef'ed away.
// Here: The 2d-version of the macros
#define UG_DIM 2
#include "uggrid/ugwrapper.hh"
#undef UG_DIM

// UG defines a whole load of preprocessor macros.  ug_undefs.hh undefines
// them all, so we don't get name clashes.
#ifdef UG_LGMDOMAIN
#include "uggrid/ug_undefs_lgm_seq.hh"
#else
#include "uggrid/ug_undefs.hh"
#endif
#undef _2

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

#define _3
#include "uggrid/ugincludes.hh"

// Wrap a few large UG macros by functions before they get undef'ed away.
// This time it's the 3d-versions.
#define UG_DIM 3
#include "uggrid/ugwrapper.hh"
#undef UG_DIM

// undef all macros defined by UG
#ifdef UG_LGMDOMAIN
#include "uggrid/ug_undefs_lgm_seq.hh"
#else
#include "uggrid/ug_undefs.hh"
#endif

#undef _3
#undef FOR_DUNE

// The components of the UGGrid interface
#include "uggrid/uggridgeometry.hh"
#include "uggrid/uggridlocalgeometry.hh"
#include "uggrid/uggridentity.hh"
#include "uggrid/uggridentitypointer.hh"
#include "uggrid/uggridentityseed.hh"
#include "uggrid/uggridintersections.hh"
#include "uggrid/uggridintersectioniterators.hh"
#include "uggrid/uggridleveliterator.hh"
#include "uggrid/uggridleafiterator.hh"
#include "uggrid/uggridhieriterator.hh"
#include "uggrid/uggridindexsets.hh"

// Not needed here, but included for user convenience
#include "uggrid/uggridfactory.hh"

#ifdef ModelP
namespace Dune {

  // converts the UG speak message buffers to DUNE speak and vice-versa
  template <class DataHandle, int GridDim, int codim>
  class UGMessageBufferBase {
  protected:
    typedef UGMessageBufferBase<DataHandle, GridDim, codim>  ThisType;
    typedef UGGrid<GridDim>                              GridType;
    typedef typename DataHandle::DataType DataType;

    enum {
      dim = GridDim
    };

    UGMessageBufferBase(void *ugData)
    {
      ugData_ = static_cast<char*>(ugData);
    };

  public:
    void write(const DataType &t)
    { this->writeRaw_<DataType>(t);  }

    void read(DataType &t)
    { this->readRaw_<DataType>(t);  }

  protected:
    friend class Dune::UGGrid<dim>;

    template <class ValueType>
    void writeRaw_(const ValueType &v)
    {
      *reinterpret_cast<ValueType*>(ugData_) = v;
      ugData_ += sizeof(ValueType);
    }

    template <class ValueType>
    void readRaw_(ValueType &v)
    {
      v = *reinterpret_cast<ValueType*>(ugData_);
      ugData_ += sizeof(ValueType);
    }

    // called by DDD_IFOneway to serialize the data structure to
    // be send
    static int ugGather_(typename UG_NS<dim>::DDD_OBJ obj, void* data)
    {
      if (codim == 0) {
        UGMakeableEntity<0, dim, UGGrid<dim> > e((typename UG_NS<dim>::Element*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Element*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.template writeRaw_<unsigned>(duneDataHandle_->size(e));
          duneDataHandle_->gather(msgBuf, e);
        }
      }
      else if (codim == dim) {
        UGMakeableEntity<dim, dim, Dune::UGGrid<dim> > e((typename UG_NS<dim>::Node*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Node*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.template writeRaw_<unsigned>(duneDataHandle_->size(e));
          duneDataHandle_->gather(msgBuf, e);
        }
      }
      else if (codim == dim - 1) {
        UGMakeableEntity<dim-1, dim, Dune::UGGrid<dim> > e((typename UG_NS<dim>::Edge*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Edge*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.template writeRaw_<unsigned>(duneDataHandle_->size(e));
          duneDataHandle_->gather(msgBuf, e);
        }
      }
      else {
        DUNE_THROW(GridError,
                   "Only node and element wise "
                   "communication is currently "
                   "supported by UGGrid");
      }

      return 0;
    }

    // called by DDD_IFOneway to deserialize the data structure
    // that has been received
    static int ugScatter_(typename UG_NS<dim>::DDD_OBJ obj, void* data)
    {
      if (codim == 0) {
        typedef UGMakeableEntity<0, dim, UGGrid<dim> > Entity;
        Entity e((typename UG_NS<dim>::Element*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Element*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          int n;
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.readRaw_(n);
          else
            n = duneDataHandle_->template size<Entity>(e);
          if (n > 0)
            duneDataHandle_->template scatter<ThisType, Entity>(msgBuf, e, n);
        }
      }
      else if (codim == dim) {
        typedef UGMakeableEntity<dim, dim, Dune::UGGrid<dim> > Entity;
        Entity e((typename UG_NS<dim>::Node*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Node*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          int n;
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.readRaw_(n);
          else
            n = duneDataHandle_->template size<Entity>(e);
          if (n > 0)
            duneDataHandle_->template scatter<ThisType, Entity>(msgBuf, e, n);
        }
      }
      else if (codim == dim - 1) {        // !!!ALEX!!! Is it possible to send codim 1 in UG<2>?
        typedef UGMakeableEntity<dim-1, dim, Dune::UGGrid<dim> > Entity;
        Entity e((typename UG_NS<dim>::Edge*)obj);
        // safety check to only communicate what is needed
        if ((level == -1 && UG_NS<dim>::isLeaf((typename UG_NS<dim>::Edge*)obj)) || e.level() == level)
        {
          ThisType msgBuf(static_cast<DataType*>(data));
          int n;
          if (!duneDataHandle_->fixedsize(dim, codim))
            msgBuf.readRaw_(n);
          else
            n = duneDataHandle_->template size<Entity>(e);
          if (n > 0)
            duneDataHandle_->template scatter<ThisType, Entity>(msgBuf, e, n);
        }
      }
      else {
        DUNE_THROW(GridError,
                   "Only node and element wise "
                   "communication is currently "
                   "supported by UGGrid");
      }

      return 0;
    }
    static DataHandle *duneDataHandle_;

    static int level;

    char              *ugData_;
  };

  template <class DataHandle, int GridDim, int codim>
  class UGMessageBuffer
    : public UGMessageBufferBase<DataHandle, GridDim, codim>
  {
    typedef typename DataHandle::DataType DataType;
    typedef UGMessageBufferBase<DataHandle, GridDim, codim> Base;
    enum { dim = GridDim };

  protected:
    friend class Dune::UGGrid<dim>;

    UGMessageBuffer(void *ugData)
      : Base(ugData)
    {}

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize_(const GridView &gv)
    {
      if (Base::duneDataHandle_->fixedsize(dim, codim)) {
        return sizeof(DataType)
               * Base::duneDataHandle_->size(*gv.template begin<codim,InteriorBorder_Partition>());
      }

      typedef typename GridView::template Codim<codim>::Entity Entity;

      // iterate over all entities, find the maximum size for
      // the current rank
      int maxSize = 0;
      typedef typename
      GridView
      ::template Codim<codim>
      ::template Partition<Dune::All_Partition>
      ::Iterator Iterator;
      Iterator it = gv.template begin<codim, Dune::All_Partition>();
      const Iterator endIt = gv.template end<codim, Dune::All_Partition>();
      for (; it != endIt; ++it) {
        maxSize = std::max((int) maxSize,
                           (int) Base::duneDataHandle_->size(*it));
      }

      // find maximum size for all ranks
      maxSize = MPIHelper::getCollectiveCommunication().max(maxSize);
      if (!maxSize)
        return 0;

      // add the size of an unsigned integer to the actual
      // buffer size. (we somewhere have to store the actual
      // number of objects for each entity.)
      return sizeof(unsigned) + sizeof(DataType)*maxSize;
    }
  };

  template <class DataHandle>
  class UGMessageBuffer<DataHandle, 2, 1>
    : public UGMessageBufferBase<DataHandle, 2, 1>
  {
    enum {codim = 1,
          GridDim = 2,
          dim = 2};
    typedef typename DataHandle::DataType DataType;
    typedef UGMessageBufferBase<DataHandle, GridDim, codim> Base;
  protected:
    friend class Dune::UGGrid<dim>;

    UGMessageBuffer(void *ugData)
      : Base(ugData)
    {}

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize_(const GridView &gv)
    {
      if (Base::duneDataHandle_->fixedsize(dim, codim)) {
        typedef typename GridView::template Codim<0>::Entity Element;
        const Element& element = *gv.template begin<0, InteriorBorder_Partition>();
        return sizeof(DataType)
               * Base::duneDataHandle_->size(element.template subEntity<codim>(0));
      }

      DUNE_THROW(GridError, "Only fixedsize implemented");
    }
  };

  template <class DataHandle>
  class UGMessageBuffer<DataHandle, 3, 2>
    : public UGMessageBufferBase<DataHandle, 3, 2>
  {
    enum {codim = 2,
          GridDim = 3,
          dim = 3};
    typedef typename DataHandle::DataType DataType;
    typedef UGMessageBufferBase<DataHandle, GridDim, codim> Base;
  protected:
    friend class Dune::UGGrid<dim>;

    UGMessageBuffer(void *ugData)
      : Base(ugData)
    {}

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize_(const GridView &gv)
    {
      if (Base::duneDataHandle_->fixedsize(dim, codim)) {
        typedef typename GridView::template Codim<0>::Entity Element;
        const Element& element = *gv.template begin<0, InteriorBorder_Partition>();
        return sizeof(DataType)
               * Base::duneDataHandle_->size(element.template subEntity<codim>(0));
      }

      DUNE_THROW(GridError, "Only fixedsize implemented");
    }
  };

}   // end namespace Dune

template <class DataHandle, int GridDim, int codim>
DataHandle *Dune::UGMessageBufferBase<DataHandle,GridDim,codim>::duneDataHandle_ = 0;

template <class DataHandle, int GridDim, int codim>
int Dune::UGMessageBufferBase<DataHandle,GridDim,codim>::level = -1;
#endif // ModelP

namespace Dune {

#ifdef ModelP
  template <int dim>
  class CollectiveCommunication<Dune::UGGrid<dim> > :
    public CollectiveCommunication< Dune::MPIHelper::MPICommunicator >
  {
    typedef CollectiveCommunication< Dune::MPIHelper::MPICommunicator > ParentType;
  public:
    CollectiveCommunication()
      : ParentType(MPIHelper::getCommunicator())
    {}
  };
#endif // ModelP

  template<int dim, int dimworld>
  struct UGGridFamily
  {
    typedef GridTraits<dim,dimworld,Dune::UGGrid<dim>,
        UGGridGeometry,
        UGGridEntity,
        UGGridEntityPointer,
        UGGridLevelIterator,
        UGGridLeafIntersection,
        UGGridLevelIntersection,
        UGGridLeafIntersectionIterator,
        UGGridLevelIntersectionIterator,
        UGGridHierarchicIterator,
        UGGridLeafIterator,
        UGGridLevelIndexSet< const UGGrid<dim> >,
        UGGridLeafIndexSet< const UGGrid<dim> >,
        UGGridIdSet< const UGGrid<dim>, false >,
        unsigned int,
        UGGridIdSet< const UGGrid<dim>, true >,
        unsigned int,
        CollectiveCommunication<Dune::UGGrid<dim> >,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
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
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief Provides the meshes of the finite element toolbox UG.
     \brief (http://atlas.gcsc.uni-frankfurt.de/~ug/).
     \ingroup GridImplementations

     This is the implementation of the grid interface
     using the UG grid management system.  UG provides conforming grids
     in two and three space dimensions.  The grids can be mixed, i.e.
     2d grids can contain triangles and quadrilaterals and 3d grids can
     contain tetrahedra and hexahedra and also pyramids and prisms.
     The grid refinement rules are very flexible.  Local adaptive red/green
     refinement is the default, but a special method in the UGGrid class
     allows you to directly access a number of anisotropic refinements rules.
     Last but not least, the UG grid manager is completely parallelized,
     and you can use boundaries parametrized by either analytical expressions
     or high-resolution piecewise linear surfaces.

     To use this module you need the UG library.  See the
     DUNE installation notes
     on how to obtain and install it.

     In your %Dune application, you can now instantiate objects of the
     type UGGrid<2> or UGGrid<3>.  You can have more than one, if
     you choose.  It is even possible to have 2d and 3d grids at the same
     time, even though the original UG system never intended to support
     this!

     See the documentation for the factory class GridFactory<UGGrid<dimworld> >
     to learn how to create UGGrid objects.

     Please send any questions, suggestions, or bug reports to the Dune mailing list
     dune@dune-project.org

     For installation instructions see http://www.dune-project.org/external_libraries/install_ug.html .
   */
  template <int dim>
  class UGGrid : public GridDefaultImplementation  <dim, dim, double, UGGridFamily<dim,dim> >
  {
    friend class UGGridGeometry<0,dim,const UGGrid<dim> >;
    friend class UGGridGeometry<dim,dim,const UGGrid<dim> >;
    friend class UGGridGeometry<1,2,const UGGrid<dim> >;
    friend class UGGridGeometry<2,3,const UGGrid<dim> >;

    friend class UGGridEntity <0,dim,const UGGrid<dim> >;
    friend class UGGridEntity <dim,dim,const UGGrid<dim> >;
    friend class UGGridHierarchicIterator<const UGGrid<dim> >;
    friend class UGGridLeafIntersection<const UGGrid<dim> >;
    friend class UGGridLevelIntersection<const UGGrid<dim> >;
    friend class UGGridLeafIntersectionIterator<const UGGrid<dim> >;
    friend class UGGridLevelIntersectionIterator<const UGGrid<dim> >;

    friend class UGGridLevelIndexSet<const UGGrid<dim> >;
    friend class UGGridLeafIndexSet<const UGGrid<dim> >;
    friend class UGGridIdSet<const UGGrid<dim>, false >;
    friend class UGGridIdSet<const UGGrid<dim>, true >;

    friend class GridFactory<UGGrid<dim> >;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLeafIterator;
    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;
    template <int codim_, class GridImp_>
    friend class UGGridEntityPointer;

    /** \brief UGGrid is only implemented for 2 and 3 dimension */
    dune_static_assert(dim==2 || dim==3, "Use UGGrid only for 2d and 3d!");

    // The different instantiations are mutual friends so they can access
    // each others numOfUGGrids field
    friend class UGGrid<2>;
    friend class UGGrid<3>;

    //**********************************************************
    // The Interface Methods
    //**********************************************************
  public:
    //! type of the used GridFamily for this grid
    typedef UGGridFamily<dim,dim>  GridFamily;

    // the Traits
    typedef typename UGGridFamily<dim,dim>::Traits Traits;

    //! The type used to store coordinates
    typedef UG::DOUBLE ctype;

    /** \brief Constructor with control over UG's memory requirements
     *
     * \param heapSize The size of UG's internal memory in megabytes for this grid.
     *
     * \note The heapsize provided here is stored as the default heapsize.
     */
    UGGrid(unsigned int heapSize) DUNE_DEPRECATED;

    /** \brief Default constructor
     *
     * Uses the default heapsize, which can be set using the static method
     * setDefaultHeapSize() (or by calling the contructor which take the
     * heapsize as an argument).
     */
    UGGrid();

    //! Destructor
    ~UGGrid();

  private:
    /** \brief Common part of the constructors
     */
    void init();

  public:
    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxlevel with 0 the coarsest level.
    int maxLevel() const;

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const;

    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const;

    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const;

    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return typename Traits::template Codim<codim>::template Partition<All_Partition>::LeafIterator(*this);
    }

    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return UGGridLeafIterator<codim,All_Partition, const UGGrid<dim> >();
    }

    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator(*this);
    }

    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return UGGridLeafIterator<codim,PiType, const UGGrid<dim> >();
    }

    /** \brief Create an EntityPointer from an EntitySeed */
    template <int codim>
    static typename Traits::template Codim<codim>::EntityPointer
    entityPointer(const UGGridEntitySeed<codim, const UGGrid<dim> >& seed)
    {
      return typename Traits::template Codim<codim>::EntityPointer(UGGridEntityPointer<codim,const UGGrid<dim> >(seed.target(),seed.gridImp()));
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
      return globalIdSet_;
    }

    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const
    {
      return localIdSet_;
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

       \param e Pointer to the element to be marked for refinement
       \param rule One of the UG refinement rules
       \param side If rule==UG::%D2::%BLUE (one quadrilateral is split into two rectangles)
       you can choose the orientation of the cut by setting side==0 or side==1
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

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return 0;
    }

    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return (codim==0) ? 1 : 0;
    }

    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return 0;
    }

    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return (codim==0) ? 1 : 0;
    }

    /** \brief Default load balancing.
        \bug The return value is always 'true'
        \return true if the grid has changed
     */
    bool loadBalance() {
      return loadBalance(0,0,2,32,1);
    }

    /** \brief Re-balances the load each process has to handle for a parallel grid,
        the DataHandle data works like the data handle for the communicate
        methods. If grid has changed , true is returned.
        \bug Not implemented yet!
     */
    template<class DataHandle>
    bool loadBalance (DataHandle& data)
    {
      DUNE_THROW(NotImplemented, "load balancing with data attached");
    }

    /** \brief Distributes this grid over the available nodes in a distributed machine

       If you want the UG default for the parameters pick
       <ul>
       <li>strategy = 0</li>
       <li>minlevel = 1</li>
       <li>depth = 2</li>
       <li>maxlevel = 32 </li>
       <li>minelement = 1</li>
       </ul>

       \bug The return value is always 'true'

       \param minlevel The coarsest grid level that gets distributed
       \param maxlevel does currently get ignored
     */
    bool loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement);

    /** \brief The communication interface for all codims on a given level
       @param dataHandle type used to gather/scatter data in and out of the message buffer
       @param iftype one of the predifined interface types, throws error if it is not implemented
       @param dir choose beetween forward and backward communication
       @param level communicate for entities on the given level

       Implements a generic communication function sending an object of type P for each entity
       in the intersection of two processors. P has two methods gather and scatter that implement
       the protocol. Therefore P is called the "protocol class".
     */
    template<class DataHandle>
    void communicate (DataHandle& dataHandle, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
#ifdef ModelP
      typedef typename UGGrid::LevelGridView LevelGridView;

      for (int curCodim = 0; curCodim <= dim; ++curCodim) {
        if (!dataHandle.contains(dim, curCodim))
          continue;

        if (curCodim == 0)
          communicateUG_<LevelGridView, DataHandle, 0>(this->levelView(level), level, dataHandle, iftype, dir);
        else if (curCodim == dim)
          communicateUG_<LevelGridView, DataHandle, dim>(this->levelView(level), level, dataHandle, iftype, dir);
        else if (curCodim == dim - 1)
          communicateUG_<LevelGridView, DataHandle, dim-1>(this->levelView(level), level, dataHandle, iftype, dir);
        else
          DUNE_THROW(NotImplemented,
                     className(*this) << "::communicate(): Only "
                     "supported for codim=0 and "
                     "codim=dim(=" << dim << "), but "
                     "codim=" << curCodim << " was requested");
      }
#endif // ModelP
    }

    /** \brief The communication interface for all codims on the leaf level
       @param dataHandle type used to gather/scatter data in and out of the message buffer
       @param iftype one of the predifined interface types, throws error if it is not implemented
       @param dir choose beetween forward and backward communication

       Implements a generic communication function sending an object of type P for each entity
       in the intersection of two processors. P has two methods gather and scatter that implement
       the protocol. Therefore P is called the "protocol class".
     */
    template<class DataHandle>
    void communicate(DataHandle& dataHandle,
                     InterfaceType iftype,
                     CommunicationDirection dir) const
    {
#ifdef ModelP
      typedef typename UGGrid::LeafGridView LeafGridView;

      for (int curCodim = 0; curCodim <= dim; ++curCodim) {
        if (!dataHandle.contains(dim, curCodim))
          continue;
        int level = -1;
        if (curCodim == 0)
          communicateUG_<LeafGridView, DataHandle, 0>(this->leafView(), level, dataHandle, iftype, dir);
        else if (curCodim == dim)
          communicateUG_<LeafGridView, DataHandle, dim>(this->leafView(), level, dataHandle, iftype, dir);
        else if (curCodim == dim - 1)   // !!!ALEX!!! Is it possible to send codim 1 in UG<2>?
        {
          communicateUG_<LeafGridView, DataHandle, dim-1>(this->leafView(), level, dataHandle, iftype, dir);
        }
        else
          DUNE_THROW(NotImplemented,
                     className(*this) << "::communicate(): Only "
                     "supported for codim=0 and "
                     "codim=dim(=" << dim << "), but "
                     "codim=" << curCodim << " was requested");
      }
#endif // ModelP
    }

    /** the collective communication */
    const CollectiveCommunication<UGGrid>& comm () const
    {
      return ccobj_;
    }

  protected:
#ifdef ModelP
    template <class GridView, class DataHandle, int codim>
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

      std::vector<typename UG_NS<dim>::DDD_IF> ugIfs;
      findDDDInterfaces_(ugIfs, iftype, codim);

      unsigned bufSize = UGMsgBuf::ugBufferSize_(gv);
      if (!bufSize)
        return;     // we don't need to communicate if we don't have any data!
      for (unsigned i=0; i < ugIfs.size(); ++i)
        UG_NS<dim>::DDD_IFOneway(ugIfs[i],
                                 ugIfDir,
                                 bufSize,
                                 &UGMsgBuf::ugGather_,
                                 &UGMsgBuf::ugScatter_);
    }

    void findDDDInterfaces_(std::vector<typename UG_NS<dim>::DDD_IF > &dddIfaces,
                            InterfaceType iftype,
                            int codim) const
    {
      dddIfaces.clear();
      switch (codim)
      {
      case 0 :
        switch (iftype) {
        case InteriorBorder_InteriorBorder_Interface :
          // do not communicate anything: Elements can not be in
          // the interior of two processes at the same time
          return;
        case InteriorBorder_All_Interface :
          // The communicated elements are in the sender's
          // interior and it does not matter what they are for
          // the receiver
          dddIfaces.push_back(UG_NS<dim>::ElementVHIF());
          return;
        case All_All_Interface :
          // It does neither matter what the communicated
          // elements are for sender nor for the receiver. If
          // they are seen by these two processes, data is send
          // and received.
          dddIfaces.push_back(UG_NS<dim>::ElementSymmVHIF());
          return;
        default :
          DUNE_THROW(GridError,
                     "Element communication not supported for "
                     "interfaces of type  "
                     << iftype);
        }

      case dim :
        switch (iftype)
        {
        case InteriorBorder_InteriorBorder_Interface :
          dddIfaces.push_back(UG_NS<dim>::BorderNodeSymmIF());
          return;
        case InteriorBorder_All_Interface :
          dddIfaces.push_back(UG_NS<dim>::BorderNodeSymmIF());
          dddIfaces.push_back(UG_NS<dim>::NodeIF());
          return;
        case All_All_Interface :
          dddIfaces.push_back(UG_NS<dim>::NodeAllIF());
          return;
        default :
          DUNE_THROW(GridError,
                     "Node communication not supported for "
                     "interfaces of type  "
                     << iftype);
        }

      case dim-1 :
        switch (iftype)
        {
        case InteriorBorder_InteriorBorder_Interface :
          dddIfaces.push_back(UG_NS<dim>::BorderEdgeSymmIF());
          return;
        case InteriorBorder_All_Interface :
          dddIfaces.push_back(UG_NS<dim>::BorderEdgeSymmIF());
          dddIfaces.push_back(UG_NS<dim>::EdgeIF());
          return;
        default :
          DUNE_THROW(GridError,
                     "Edge communication not supported for "
                     "interfaces of type  "
                     << iftype);
        }

      default :
        DUNE_THROW(GridError,
                   "Communication for codim "
                   << codim
                   << " entities is not yet supported "
                   << " by the DUNE UGGrid interface!");
      }
    };
#endif // ModelP

  public:
    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    /** @name Coarse Grid Creation Methods */
    /*@{*/

    /** \brief When UGGrid has been configured to use the LGM domain manager,
        this routine sets up a grid from an LGM and an NG file
     */
    void createLGMGrid(const std::string& name);

    /*@}*/

    /** \brief Rudimentary substitute for a hierarchic iterator on faces
        \param e, elementSide Grid face specified by an element and one of its sides
        \param maxl The finest level that should be traversed by the iterator
        \param[out] childElements For each subface: element index, elementSide, and level
     */
    void getChildrenOfSubface(const typename Traits::template Codim<0>::EntityPointer & e,
                              int elementSide,
                              int maxl,
                              std::vector<typename Traits::template Codim<0>::EntityPointer>& childElements,
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

    /** \brief Sets the default heap size
     *
     * UGGrid keeps an internal heap to allocate memory from, which must be
     * specified on grid creation (at the latest).  This sets the default heap
     * size, which is used when no heap size is given to the constructor.
     */
    static void setDefaultHeapSize(unsigned size) {
      heapSize_ = size;
    }

    /** \brief Sets a vertex to a new position

       Changing a vertex' position changes its position on all grid levels!*/
    void setPosition(const typename Traits::template Codim<dim>::EntityPointer& e,
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

    /** \brief The collective communication object. */
    CollectiveCommunication<UGGrid> ccobj_;

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
    std::vector<UGGridLevelIndexSet<const UGGrid<dim> >*> levelIndexSets_;

    UGGridLeafIndexSet<const UGGrid<dim> > leafIndexSet_;

    UGGridIdSet<const UGGrid<dim>, false > globalIdSet_;

    UGGridIdSet<const UGGrid<dim>, true > localIdSet_;

    //! The type of grid refinement currently in use
    RefinementType refinementType_;

    //! The type of grid refinement closure currently in use
    ClosureType closureType_;

    /** \brief Number of UGGrids currently in use.
     *
     * This counts the number of UGGrids currently instantiated.  All
     * constructors of UGGrid look at this variable.  If it zero, they
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

    /** \brief The size of UG's internal heap in megabytes
     *
     * It is handed over to UG for each new multigrid.
     */
    static unsigned int heapSize_;

    /** \brief The classes implementing the geometry of the boundary segments, if requested */
    std::vector<shared_ptr<BoundarySegment<dim> > > boundarySegments_;

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

    /** \brief UGGrid has codim=0 entities (elements)
       \ingroup UGGrid
     */
    template<int dim>
    struct hasEntity< UGGrid<dim>, 0>
    {
      static const bool v = true;
    };

    /** \brief UGGrid has codim=dim entities (vertices)
       \ingroup UGGrid
     */
    template<int dim>
    struct hasEntity< UGGrid<dim>, dim>
    {
      static const bool v = true;
    };

    /** \brief UGGrid is parallel
       \ingroup UGGrid
     */
    template<int dim>
    struct isParallel< UGGrid<dim> >
    {
#ifdef ModelP
      static const bool v = true;
#else // !ModelP
      static const bool v = false;
#endif // !ModelP
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

#endif   // HAVE_UG
#endif   // DUNE_UGGRID_HH
