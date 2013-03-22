// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_YASPGRID_HH
#define DUNE_YASPGRID_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>

// either include stdint.h or provide fallback for uint8_t
#if HAVE_STDINT_H
#include <stdint.h>
#else
typedef unsigned char uint8_t;
#endif

#include <dune/grid/common/grid.hh>     // the grid base classes
#include <dune/grid/yaspgrid/grids.hh>  // the yaspgrid base classes
#include <dune/grid/common/capabilities.hh> // the capabilities
#include <dune/common/misc.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/datahandleif.hh>


#if HAVE_MPI
#include <dune/common/parallel/mpicollectivecommunication.hh>
#endif

/*! \file yaspgrid.hh
   YaspGrid stands for yet another structured parallel grid.
   It will implement the dune grid interface for structured grids with codim 0
   and dim, with arbitrary overlap, parallel features with two overlap
   models, periodic boundaries and fast a implementation allowing on-the-fly computations.
 */

namespace Dune {

  //************************************************************************
  /*! define name for floating point type used for coordinates in yaspgrid.
     You can change the type for coordinates by changing this single typedef.
   */
  typedef double yaspgrid_ctype;

  /* some sizes for building global ids
   */
  const int yaspgrid_dim_bits = 24; // bits for encoding each dimension
  const int yaspgrid_level_bits = 6; // bits for encoding level number
  const int yaspgrid_codim_bits = 4; // bits for encoding codimension


  //************************************************************************
  // forward declaration of templates

  template<int dim>                             class YaspGrid;
  template<int mydim, int cdim, class GridImp>  class YaspGeometry;
  template<int codim, int dim, class GridImp>   class YaspEntity;
  template<int codim, class GridImp>            class YaspEntityPointer;
  template<int codim, class GridImp>            class YaspEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp> class YaspLevelIterator;
  template<class GridImp>            class YaspIntersectionIterator;
  template<class GridImp>            class YaspIntersection;
  template<class GridImp>            class YaspHierarchicIterator;
  template<class GridImp>            class YaspIndexSet;
  template<class GridImp>            class YaspGlobalIdSet;

  namespace FacadeOptions
  {

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, const YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

  }

} // namespace Dune

#include <dune/grid/yaspgrid/yaspgridgeometry.hh>
#include <dune/grid/yaspgrid/yaspgridentity.hh>
#include <dune/grid/yaspgrid/yaspgridintersection.hh>
#include <dune/grid/yaspgrid/yaspgridintersectioniterator.hh>
#include <dune/grid/yaspgrid/yaspgridhierarchiciterator.hh>
#include <dune/grid/yaspgrid/yaspgridentityseed.hh>
#include <dune/grid/yaspgrid/yaspgridentitypointer.hh>
#include <dune/grid/yaspgrid/yaspgridleveliterator.hh>
#include <dune/grid/yaspgrid/yaspgridindexsets.hh>
#include <dune/grid/yaspgrid/yaspgrididset.hh>

namespace Dune {

  template<int dim>
  struct YaspGridFamily
  {
#if HAVE_MPI
    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
    typedef CollectiveCommunication<Dune::YaspGrid<dim> > CCType;
#endif

    typedef GridTraits<dim,                                     // dimension of the grid
        dim,                                                    // dimension of the world space
        Dune::YaspGrid<dim>,
        YaspGeometry,YaspEntity,
        YaspEntityPointer,
        YaspLevelIterator,                                      // type used for the level iterator
        YaspIntersection,              // leaf  intersection
        YaspIntersection,              // level intersection
        YaspIntersectionIterator,              // leaf  intersection iter
        YaspIntersectionIterator,              // level intersection iter
        YaspHierarchicIterator,
        YaspLevelIterator,                                      // type used for the leaf(!) iterator
        YaspIndexSet< const YaspGrid< dim > >,                  // level index set
        YaspIndexSet< const YaspGrid< dim > >,                  // leaf index set
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        CCType,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        YaspEntitySeed>
    Traits;
  };

  template<int dim, int codim>
  struct YaspCommunicateMeta {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,codim))
      {
        DUNE_THROW(GridError, "interface communication not implemented");
      }
      YaspCommunicateMeta<dim,codim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,dim> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,dim))
        g.template communicateCodim<DataHandle,dim>(data,iftype,dir,level);
      YaspCommunicateMeta<dim,dim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,0> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,0))
        g.template communicateCodim<DataHandle,0>(data,iftype,dir,level);
    }
  };


  //************************************************************************
  /*!
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief Provides a distributed structured cube mesh.
     \ingroup GridImplementations

     YaspGrid stands for yet another structured parallel grid.
     It implements the dune grid interface for structured grids with codim 0
     and dim, with arbitrary overlap (including zero),
     periodic boundaries and fast implementation allowing on-the-fly computations.

     \tparam dim The dimension of the grid and its surrounding world

     \par History:
     \li started on July 31, 2004 by PB based on abstractions developed in summer 2003
   */
  template<int dim>
  class YaspGrid :
    public GridDefaultImplementation<dim,dim,yaspgrid_ctype,YaspGridFamily<dim> >,
    public MultiYGrid<dim,yaspgrid_ctype>
  {
    typedef const YaspGrid<dim> GridImp;

    void init()
    {
      setsizes();
      indexsets.push_back( make_shared< YaspIndexSet<const YaspGrid<dim> > >(*this,0) );
      boundarysegmentssize();
    }

    void boundarysegmentssize()
    {
      // sizes of local macro grid
      const FieldVector<int, dim> & size = MultiYGrid<dim,ctype>::begin().cell_overlap().size();
      FieldVector<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((MultiYGrid<dim,ctype>::begin().cell_overlap().origin(i)
              == MultiYGrid<dim,ctype>::begin().cell_global().origin(i))+
             (MultiYGrid<dim,ctype>::begin().cell_overlap().origin(i) +
                    MultiYGrid<dim,ctype>::begin().cell_overlap().size(i)
                    == MultiYGrid<dim,ctype>::begin().cell_global().origin(i) +
                    MultiYGrid<dim,ctype>::begin().cell_global().size(i)));
        }
      }
      nBSegments = 0;
      for (int k=0; k<dim; k++)
      {
        int offset = 1;
        for (int l=0; l<dim; l++)
        {
          if (l==k) continue;
          offset *= size[l];
        }
        nBSegments += sides[k]*offset;
      }
    }

  public:

    using MultiYGrid<dim,yaspgrid_ctype>::defaultLoadbalancer;

    //! define type used for coordinates in grid module
    typedef yaspgrid_ctype ctype;

    // define the persistent index type
    typedef bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits> PersistentIndexType;

    //! the GridFamily of this grid
    typedef YaspGridFamily<dim> GridFamily;
    // the Traits
    typedef typename YaspGridFamily<dim>::Traits Traits;

    // need for friend declarations in entity
    typedef YaspIndexSet<YaspGrid<dim> > LevelIndexSetType;
    typedef YaspIndexSet<YaspGrid<dim> > LeafIndexSetType;
    typedef YaspGlobalIdSet<YaspGrid<dim> > GlobalIdSetType;

    //! maximum number of levels allowed
    enum { MAXL=64 };

    //! shorthand for base class data types
    typedef MultiYGrid<dim,ctype> YMG;
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef typename MultiYGrid<dim,ctype>::Intersection IS;
    typedef typename std::deque<IS>::const_iterator ISIT;

    /*! Constructor for a YaspGrid, they are all forwarded to the base class
       @param comm MPI communicator where this mesh is distributed to
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::MPIHelper::MPICommunicator comm,
              Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
#if HAVE_MPI
      : YMG(comm,L,s,std::bitset<dim>(),overlap,lb), ccobj(comm),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,std::bitset<dim>(),overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      // hack: copy input bitfield (in FieldVector<bool>) into std::bitset
      for (size_t i=0; i<dim; i++)
        this->_periodic[i] = periodic[i];
      init();
    }


    /*! Constructor for a sequential YaspGrid, they are all forwarded to the base class.

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = YMG::defaultLoadbalancer())
#if HAVE_MPI
      : YMG(MPI_COMM_SELF,L,s,std::bitset<dim>(),overlap,lb), ccobj(MPI_COMM_SELF),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,std::bitset<dim>(),overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      // hack: copy input bitfield (in FieldVector<bool>) into std::bitset
      for (size_t i=0; i<dim; i++)
        this->_periodic[i] = periodic[i];
      init();
    }

    /*! Constructor for a YaspGrid, they are all forwarded to the base class
       @param comm MPI communicator where this mesh is distributed to
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::MPIHelper::MPICommunicator comm,
              Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> s,
              std::bitset<dim> periodic,
              int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
#if HAVE_MPI
      : YMG(comm,L,s,periodic,overlap,lb), ccobj(comm),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,periodic,overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      init();
    }


    /*! Constructor for a sequential YaspGrid, they are all forwarded to the base class.

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> s,
              std::bitset<dim> periodic,
              int overlap,
              const YLoadBalance<dim>* lb = YMG::defaultLoadbalancer())
#if HAVE_MPI
      : YMG(MPI_COMM_SELF,L,s,periodic,overlap,lb), ccobj(MPI_COMM_SELF),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,periodic,overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      init();
    }

    /*! Constructor for a sequential YaspGrid without periodicity

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain (lower left is always (0,...,0)
       @param elements number of cells on coarse mesh in each direction
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::array<int, dim> elements)
#if HAVE_MPI
      : YMG(MPI_COMM_SELF,L,elements,std::bitset<dim>(0),0), ccobj(MPI_COMM_SELF),
#else
      : YMG(L,elements,std::bitset<dim>(0),0),
#endif
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
    {
      init();
    }

  private:
    // do not copy this class
    YaspGrid(const YaspGrid&);

  public:

    /*! Return maximum level defined in this grid. Levels are numbered
          0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {return MultiYGrid<dim,ctype>::maxlevel();} // delegate

    //! refine the grid refCount times. What about overlap?
    void globalRefine (int refCount)
    {
      if (refCount < -maxLevel())
        DUNE_THROW(GridError, "Only " << maxLevel() << " levels left. " <<
                   "Coarsening " << -refCount << " levels requested!");
      for (int k=refCount; k<0; k++)
      {
        MultiYGrid<dim,ctype>::coarsen();
        setsizes();
        indexsets.pop_back();
      }
      for (int k=0; k<refCount; k++)
      {
        MultiYGrid<dim,ctype>::refine(keep_ovlp);
        setsizes();
        indexsets.push_back( make_shared<YaspIndexSet<const YaspGrid<dim> > >(*this,maxLevel()) );
      }
    }

    /**
       \brief set options for refinement
       @param keepPhysicalOverlap [true] keep the physical size of the overlap, [false] keep the number of cells in the overlap.  Default is [true].
     */
    void refineOptions (bool keepPhysicalOverlap)
    {
      keep_ovlp = keepPhysicalOverlap;
    }

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity to Entity that should be refined

       \return true if Entity was marked, false otherwise.

       \note
          -  On yaspgrid marking one element will mark all other elements of the level as well
          -  If refCount is lower than refCount of a previous mark-call, nothing is changed
     */
    bool mark( int refCount, const typename Traits::template Codim<0>::Entity & e )
    {
      assert(adaptActive == false);
      if (e.level() != maxLevel()) return false;
      adaptRefCount = std::max(adaptRefCount, refCount);
      return true;
    }

    /** \brief returns adaptation mark for given entity

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark, here the default value 0 is returned
     */
    int getMark ( const typename Traits::template Codim<0>::Entity &e ) const
    {
      return ( e.level() == maxLevel() ) ? adaptRefCount : 0;
    }

    //! map adapt to global refine
    bool adapt ()
    {
      globalRefine(adaptRefCount);
      return (adaptRefCount > 0);
    }

    //! returns true, if the grid will be coarsened
    bool preAdapt ()
    {
      adaptActive = true;
      adaptRefCount = comm().max(adaptRefCount);
      return (adaptRefCount < 0);
    }

    //! clean up some markers
    void postAdapt()
    {
      adaptActive = false;
      adaptRefCount = 0;
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,pitype>(level);
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return levelend<cd,pitype>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,All_Partition>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      return levelend<cd,All_Partition>(level);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
      return levelend<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,All_Partition>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      return levelend<cd,All_Partition>(maxLevel());
    }

    // \brief obtain EntityPointer from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::EntityPointer
    entityPointer(const Seed& seed) const
    {
      const int codim = Seed::codimension;
      YGLI g = MultiYGrid<dim,ctype>::begin(this->getRealImplementation(seed).level());
      switch (codim)
      {
      case 0 :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g.cell_overlap(), this->getRealImplementation(seed).coord()));
      case dim :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g.vertex_overlap(), this->getRealImplementation(seed).coord()));
      default :
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim) const
    {
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      return g.overlap();
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int codim) const
    {
      YGLI g = MultiYGrid<dim,ctype>::begin(maxLevel());
      return g.overlap();
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int codim) const
    {
      return 0;
    }

    //! number of entities per level and codim in this process
    int size (int level, int codim) const
    {
      return sizes[level][codim];
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return sizes[maxLevel()][codim];
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? sizes[level][dim-type.dim()] : 0;
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return size(maxLevel(),type);
    }

    //! \brief returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const
    {
      return nBSegments;
    }

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,level);
    }

    /*! The new communication interface

       communicate objects for all codims on the leaf grid
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,this->maxLevel());
    }

    /*! The new communication interface

       communicate objects for one codim
     */
    template<class DataHandle, int codim>
    void communicateCodim (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      // check input
      if (!data.contains(dim,codim)) return; // should have been checked outside

      // data types
      typedef typename DataHandle::DataType DataType;

      // access to grid level
      YGLI g = MultiYGrid<dim,ctype>::begin(level);

      // find send/recv lists or throw error
      const std::deque<IS>* sendlist=0;
      const std::deque<IS>* recvlist=0;
      if (codim==0) // the elements
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
          return; // there is nothing to do in this case
        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g.send_cell_interior_overlap();
          recvlist = &g.recv_cell_overlap_interior();
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface || iftype==All_All_Interface)
        {
          sendlist = &g.send_cell_overlap_overlap();
          recvlist = &g.recv_cell_overlap_overlap();
        }
      }
      if (codim==dim) // the vertices
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
        {
          sendlist = &g.send_vertex_interiorborder_interiorborder();
          recvlist = &g.recv_vertex_interiorborder_interiorborder();
        }

        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g.send_vertex_interiorborder_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_interiorborder();
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface)
        {
          sendlist = &g.send_vertex_overlap_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_overlap();
        }
        if (iftype==All_All_Interface)
        {
          sendlist = &g.send_vertex_overlapfront_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_overlapfront();
        }
      }

      // change communication direction?
      if (dir==BackwardCommunication)
        std::swap(sendlist,recvlist);

      int cnt;

      // Size computation (requires communication if variable size)
      std::vector<int> send_size(sendlist->size(),-1);    // map rank to total number of objects (of type DataType) to be sent
      std::vector<int> recv_size(recvlist->size(),-1);    // map rank to total number of objects (of type DataType) to be recvd
      std::vector<size_t*> send_sizes(sendlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be sent
      std::vector<size_t*> recv_sizes(recvlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be recvd
      if (data.fixedsize(dim,codim))
      {
        // fixed size: just take a dummy entity, size can be computed without communication
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          send_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          recv_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
      }
      else
      {
        // variable size case: sender side determines the size
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          // allocate send buffer for sizes per entitiy
          size_t *buf = new size_t[is->grid.totalsize()];
          send_sizes[cnt] = buf;

          // loop over entities and ask for size
          int i=0; size_t n=0;
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
          {
            buf[i] = data.size(*it);
            n += buf[i];
            i++;
          }

          // now we know the size for this rank
          send_size[cnt] = n;

          // hand over send request to torus class
          MultiYGrid<dim,ctype>::torus().send(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // allocate recv buffers for sizes and store receive request
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // allocate recv buffer
          size_t *buf = new size_t[is->grid.totalsize()];
          recv_sizes[cnt] = buf;

          // hand over recv request to torus class
          MultiYGrid<dim,ctype>::torus().recv(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // exchange all size buffers now
        MultiYGrid<dim,ctype>::torus().exchange();

        // release send size buffers
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          delete[] send_sizes[cnt];
          send_sizes[cnt] = 0;
          cnt++;
        }

        // process receive size buffers
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // get recv buffer
          size_t *buf = recv_sizes[cnt];

          // compute total size
          size_t n=0;
          for (int i=0; i<is->grid.totalsize(); ++i)
            n += buf[i];

          // ... and store it
          recv_size[cnt] = n;
          ++cnt;
        }
      }


      // allocate & fill the send buffers & store send request
      std::vector<DataType*> sends(sendlist->size(), static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        //      std::cout << "[" << this->comm().rank() << "] "
        //                << " send " << " dest=" << is->rank
        //                << " size=" << send_size[cnt]
        //                << std::endl;

        // allocate send buffer
        DataType *buf = new DataType[send_size[cnt]];

        // remember send buffer
        sends[cnt] = buf;

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // fill send buffer; iterate over cells in intersection
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
        for ( ; it!=tsubend; ++it)
          data.gather(mb,*it);

        // hand over send request to torus class
        MultiYGrid<dim,ctype>::torus().send(is->rank,buf,send_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // allocate recv buffers and store receive request
      std::vector<DataType*> recvs(recvlist->size(),static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        //      std::cout << "[" << this->comm().rank() << "] "
        //                << " recv " << "  src=" << is->rank
        //                << " size=" << recv_size[cnt]
        //                << std::endl;

        // allocate recv buffer
        DataType *buf = new DataType[recv_size[cnt]];

        // remember recv buffer
        recvs[cnt] = buf;

        // hand over recv request to torus class
        MultiYGrid<dim,ctype>::torus().recv(is->rank,buf,recv_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // exchange all buffers now
      MultiYGrid<dim,ctype>::torus().exchange();

      // release send buffers
      cnt=0;
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        delete[] sends[cnt];
        sends[cnt] = 0;
        cnt++;
      }

      // process receive buffers and delete them
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        // get recv buffer
        DataType *buf = recvs[cnt];

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // copy data from receive buffer; iterate over cells in intersection
        if (data.fixedsize(dim,codim))
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          size_t n=data.size(*it);
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
            data.scatter(mb,*it,n);
        }
        else
        {
          int i=0;
          size_t *sbuf = recv_sizes[cnt];
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
            data.scatter(mb,*it,sbuf[i++]);
          delete[] sbuf;
        }

        // delete buffer
        delete[] buf; // hier krachts !
        cnt++;
      }
    }

    // The new index sets from DDM 11.07.2005
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LocalIdSet& localIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      return *(indexsets[level]);
    }

    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return *indexsets.back();
    }

#if HAVE_MPI
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<MPI_Comm>& comm () const
    {
      return ccobj;
    }
#else
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<YaspGrid>& comm () const
    {
      return ccobj;
    }
#endif

  private:

#if HAVE_MPI
    CollectiveCommunication<MPI_Comm> ccobj;
#else
    CollectiveCommunication<YaspGrid> ccobj;
#endif

    std::vector< shared_ptr< YaspIndexSet<const YaspGrid<dim> > > > indexsets;
    YaspGlobalIdSet<const YaspGrid<dim> > theglobalidset;

    // number of boundary segments of the level 0 grid
    int nBSegments;

    // Index classes need access to the real entity
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspGlobalIdSet<const Dune::YaspGrid<dim> >;

    friend class Dune::YaspIntersectionIterator<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspIntersection<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspEntity<0, dim, const Dune::YaspGrid<dim> >;

    template <int codim_, class GridImp_>
    friend class Dune::YaspEntityPointer;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    template<class DT>
    class MessageBuffer {
    public:
      // Constructor
      MessageBuffer (DT *p)
      {
        a=p;
        i=0;
        j=0;
      }

      // write data to message buffer, acts like a stream !
      template<class Y>
      void write (const Y& data)
      {
        dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
        a[i++] = data;
      }

      // read data from message buffer, acts like a stream !
      template<class Y>
      void read (Y& data) const
      {
        dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
        data = a[j++];
      }

    private:
      DT *a;
      int i;
      mutable int j;
    };

    void setsizes ()
    {
      for (YGLI g=MultiYGrid<dim,ctype>::begin(); g!=MultiYGrid<dim,ctype>::end(); ++g)
      {
        // codim 0 (elements)
        sizes[g.level()][0] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g.level()][0] *= g.cell_overlap().size(i);

        // codim 1 (faces)
        if (dim>1)
        {
          sizes[g.level()][1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g.cell_overlap().size(i)+1;
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g.cell_overlap().size(j);
            sizes[g.level()][1] += s;
          }
        }

        // codim dim-1 (edges)
        if (dim>2)
        {
          sizes[g.level()][dim-1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g.cell_overlap().size(i);
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g.cell_overlap().size(j)+1;
            sizes[g.level()][dim-1] += s;
          }
        }

        // codim dim (vertices)
        sizes[g.level()][dim] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g.level()][dim] *= g.vertex_overlapfront().size(i);
      }
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelbegin (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (pitype==Ghost_Partition)
        return levelend <cd, pitype> (level);
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_interior().tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_overlap().tsubbegin());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interior().tsubbegin());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interiorborder().tsubbegin());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlap().tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlapfront().tsubbegin());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelend (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_interior().tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_overlap().tsubend());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interior().tsubend());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interiorborder().tsubend());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlap().tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlapfront().tsubend());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    int sizes[MAXL][dim+1]; // total number of entities per level and codim
    bool keep_ovlp;
    int adaptRefCount;
    bool adaptActive;
  };

  namespace Capabilities
  {

    /** \struct hasEntity
       \ingroup YaspGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup YaspGrid
     */

    /** \brief YaspGrid has only one geometry type for codim 0 entities
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasSingleGeometryType< YaspGrid<dim> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< dim > :: type :: id ;
    };

    /** \brief YaspGrid is a Cartesian grid
        \ingroup YaspGrid
     */
    template<int dim>
    struct isCartesian< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=0 entities (elements)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, 0 >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=dim entities (vertices)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, dim >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, 0 >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, dim >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is parallel
       \ingroup YaspGrid
     */
    template<int dim>
    struct isParallel< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is levelwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLevelwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is leafwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLeafwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

  }

} // end namespace


#endif
