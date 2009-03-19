// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_HH
#define DUNE_ONE_D_GRID_HH

#include <vector>
#include <list>

#include <dune/common/misc.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/grid.hh>


/** \file
 * \brief The OneDGrid class
 */

namespace Dune
{

  // forward declarations
  template<int codim, int dim, class GridImp> class OneDGridEntity;
  template<int codim, class GridImp> class OneDGridEntityPointer;
  template<int codim, PartitionIteratorType pitype, class GridImp> class OneDGridLevelIterator;

  template<int mydim, int coorddim, class GridImp>            class OneDGridGeometry;
  template<class GridImp>            class OneDGridHierarchicIterator;
  template<class GridImp, bool LeafIterator> class OneDGridIntersectionIterator;
  class OneDGrid;

  template<int codim>                        class OneDGridLevelIteratorFactory;

}

#include "onedgrid/onedgridlist.hh"
#include "onedgrid/nulliteratorfactory.hh"
#include "onedgrid/onedgridentity.hh"
#include "onedgrid/onedgridentitypointer.hh"
#include "onedgrid/onedgridgeometry.hh"
#include "onedgrid/onedintersectionit.hh"
#include "onedgrid/onedgridleveliterator.hh"
#include "onedgrid/onedgridleafiterator.hh"
#include "onedgrid/onedgridhieriterator.hh"
#include "onedgrid/onedgridindexsets.hh"

namespace Dune {

  template<int dim, int dimw>
  struct OneDGridFamily
  {
    typedef GridTraits<dim,dimw,Dune::OneDGrid,
        OneDGridGeometry,
        OneDGridEntity,
        OneDGridEntityPointer,
        OneDGridLevelIterator,
        OneDGridLeafIntersectionIterator,                // leaf  intersection
        OneDGridLevelIntersectionIterator,                // level intersection
        OneDGridLeafIntersectionIterator,                // leaf  intersection iter
        OneDGridLevelIntersectionIterator,                // level intersection iter
        OneDGridHierarchicIterator,
        OneDGridLeafIterator,
        OneDGridLevelIndexSet<const OneDGrid>,
        OneDGridLeafIndexSet<const OneDGrid>,
        OneDGridIdSet<const OneDGrid>,
        unsigned int,
        OneDGridIdSet<const OneDGrid>,
        unsigned int,
        CollectiveCommunication<Dune::OneDGrid> >
    Traits;
  };

  //**********************************************************************
  //
  // --OneDGrid
  //
  //**********************************************************************

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     Onedimensional adaptive grid
     \ingroup GridImplementations

     This implementation of the grid interface provides one-dimensional
     grids only.  No matter what the values of dim and dimworld may be,
     you'll always get a 1D-grid in a 1D-world.  Unlike SGrid, however,
     which can also be instantiated in 1D, the OneDGrid is nonuniform
     and provides local mesh refinement and coarsening.
   */
  class OneDGrid : public GridDefaultImplementation <1, 1,double,OneDGridFamily<1,1> >
  {
    // Grid and world dimension are hardwired in this grid
    enum {dim = 1};
    enum {dimworld = 1};

    friend class OneDGridLevelIteratorFactory <0>;
    friend class OneDGridLevelIteratorFactory <1>;
    template <int codim_, int dim_, class GridImp_>
    friend class OneDGridEntity;
    friend class OneDGridHierarchicIterator<OneDGrid>;
    friend class OneDGridLeafIntersectionIterator<const OneDGrid>;
    friend class OneDGridLevelIntersectionIterator<const OneDGrid>;

    friend class OneDGridLevelIndexSet<const OneDGrid>;
    friend class OneDGridLeafIndexSet<const OneDGrid>;
    friend class OneDGridIdSet<const OneDGrid>;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class OneDGridLeafIterator;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    // **********************************************************
    // The Interface Methods
    // **********************************************************

  public:

    /** \brief The type used to store coordinates

       If you ever want OneDGrid to use a different type for coordinates,
       you need to change this type and the third template argument of
       the base class.
     */
    typedef double ctype;

    /** \brief GridFamily of OneDGrid */
    typedef OneDGridFamily<dim,dimworld> GridFamily;

    //Provides the standard grid types
    typedef OneDGridFamily<dim,dimworld>::Traits Traits;

    /** \brief Constructor with an explicit set of coordinates */
    OneDGrid(const std::vector<ctype>& coords);

    /** \brief Constructor for a uniform grid */
    OneDGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary);

    //! Destructor
    ~OneDGrid();

    /** \brief Return maximum level defined in this grid.

       Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {return vertices.size()-1;}

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

    //! Iterator to first entity of given codim on leaf level
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin () const;

    //! one past the end on leaf level
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend () const;

    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const;

    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const;

    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      if (codim<0 || codim>1)
        DUNE_THROW(GridError, "There are no codim " << codim << " entities in a OneDGrid!");

      if (codim==0)
        return elements[level].size();

      return vertices[level].size();
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

    /** \brief The processor overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int overlapSize(int codim) const {
      return 0;
    }

    /** \brief The processor ghost overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int ghostSize(int codim) const {
      return 0;
    }

    /** \brief The processor overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int overlapSize(int level, int codim) const {
      return 0;
    }

    /** \brief The processor ghost overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int ghostSize(int level, int codim) const {
      return 0;
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
     * \param e EntityPointer to the entity you want to mark
     *
     * \return True, if marking was successfull
     */
    bool mark(int refCount, const Traits::Codim<0>::EntityPointer& e ) DUNE_DEPRECATED;

    /** \brief Mark entity for refinement
     *
     * \param refCount if >0 mark for refinement, if <0 mark for coarsening
     * \param e Entity to the entity you want to mark
     *
     * \return True, if marking was successfull
     */
    bool mark(int refCount, const Traits::Codim<0>::Entity& e );

    /** \brief return current adaptation marker of given entity

        \param e Entity to the entity you want to mark

        \return int current adaptation marker of entity pointer e
     */
    int getMark(const Traits::Codim<0>::EntityPointer & e ) const DUNE_DEPRECATED;

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

    /** \brief grid identification */
    std::string name () const { return "OneDGrid"; }

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

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {}

    const CollectiveCommunication &comm () const
    {
      return ccobj;
    }


  private:

    CollectiveCommunication ccobj;

    /** \brief Update all indices and ids */
    void setIndices();

    unsigned int getNextFreeId(int codim) {
      return (codim==0) ? freeElementIdCounter_++ : freeVertexIdCounter_++;
    }

    //! The type of grid refinement currently in use
    RefinementType refinementType_;

    OneDGridList<OneDEntityImp<0> >::iterator getLeftUpperVertex(const OneDEntityImp<1>* eIt);

    OneDGridList<OneDEntityImp<0> >::iterator getRightUpperVertex(const OneDEntityImp<1>* eIt);

    /** \brief Returns an iterator to the first element on the left of
        the input element which has sons.
     */
    OneDGridList<OneDEntityImp<1> >::iterator getLeftNeighborWithSon(OneDGridList<OneDEntityImp<1> >::iterator eIt);

    // The vertices of the grid hierarchy
    std::vector<OneDGridList<OneDEntityImp<0> > > vertices;
    //std::vector<std::list<OneDEntityImp<0> > > vertices;

    // The elements of the grid hierarchy
    std::vector<OneDGridList<OneDEntityImp<1> > > elements;

    // Our set of level indices
    mutable std::vector<OneDGridLevelIndexSet<const OneDGrid>* > levelIndexSets_;

    OneDGridLeafIndexSet<const OneDGrid> leafIndexSet_;

    OneDGridIdSet<const OneDGrid> idSet_;

    unsigned int freeVertexIdCounter_;

    unsigned int freeElementIdCounter_;

  }; // end Class OneDGrid

  namespace Capabilities
  {
    /** \struct hasBackupRestoreFacilities
       \ingroup OneDGrid
     */

    /** \struct IsUnstructured
       \ingroup OneDGrid
     */

    /** \brief OneDGrid has entities for all codimension
       \ingroup OneDGrid
     */
    template<int cdim>
    struct hasEntity< OneDGrid, cdim >
    {
      static const bool v = true;
    };

    /** \brief OneDGrid is not parallel
       \ingroup OneDGrid
     */
    template<>
    struct isParallel< OneDGrid >
    {
      static const bool v = false;
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

    /** \brief OneDGrid does not support hanging nodes
       \ingroup OneDGrid
     */
    template<>
    struct hasHangingNodes< OneDGrid >
    {
      static const bool v = false;
    };

  }

} // namespace Dune

#endif
