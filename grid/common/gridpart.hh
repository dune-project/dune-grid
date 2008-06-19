// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDPART_HH
#define DUNE_GRIDPART_HH

#warning "GridParts are deprecated. Use GridViews instead."
#warning "For more information see dune/grid/common/gridview.hh."

//- System includes

//- Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/deprecated.hh>

/** @file
   @author Robert Kloefkorn
   @brief Provides views of grid via grid parts, heavily used in the
   dune-fem module.
 */
/*! @addtogroup GridPart
    @ingroup Grid

   @section GridPart1 What is a GridPart ?
   <!--============================-->
   A Dune::GridPart define a view for a given grid, which can basically be seen
   as a container for entities of different co-dimensions. The most
   prominent example of such views are Leafgrids and Levelgrids - but
   also a set of entities satisfying a certain constraint, e.g., belonging
   to a certain domain, can be accessed using the GridParts interface.

   The Dune::Grid instance is passed to the GridPart which then
   provides iterators over the entities in the view; also intersection iterators
   suitable for element in the view can be obtained from the GridPart.
   Finally the GridPart provides a Dune::IndexSet with indices for all
   entities in the view. For parallel computations the GridPart interface
   also provides the suitable communication method.

   @section GridPart3 GridPart Interface and available Implementations
   <!--==========================-->

   This interface is implemented by the class template Dune::GridPartInterface.
   For a full documentation see the description of this class.
   A short list of the most important methods is:
   - The pair of methods begin()/end() provide begin and end iterators
    for the entities.
   - The methods ibegin(entity)/iend(entity)
    begin and end intersection iterators for a given entity.
   - The index set is accessed via the method indexSet() and the underlying
    grid using the method grid().
   .

   For a level or leaf view of the grid use the implementations
   Dune::LevelGridPart and Dune::LeafGridPart, respectively. A first
   implementation of a view using a general filter is available
   in the dune-fem package (Dune::FilteredGridPart).
 */

namespace Dune {
  /**
   * @addtogroup GridPart
   *
   * @{
   */

  // Forward declarations
  template <class GridImp, PartitionIteratorType pitype>
  class LevelGridPartTraits;
  template <class GridImp, PartitionIteratorType pitype>
  class LeafGridPartTraits;
  template <class GridImp, PartitionIteratorType pitype>
  class HierarchicGridPartTraits;

  //! \brief Interface for the GridPart classes
  //! A GridPart class allows to access only a specific subset of a grid's
  //! entities. A GridPart implementation provides the corresponding index set
  //! and a begin/end iterator pair for accessing those entities, the
  //! corresponding intersection iterators and a appropriate communication
  //! method.
  //! GridParts are used to parametrize spaces (see DiscreteFunctionSpaceDefault [in dune-fem]).
  template <class GridPartTraits>
  class GridPartInterface {
  public:
    //! \brief Type of the Traits
    typedef GridPartTraits Traits;

    //! \brief Type of the implementation
    typedef typename GridPartTraits::GridPartType GridPartType;


    //! \brief type of Grid implementation
    typedef typename GridPartTraits::GridType GridType;

    //! \brief Index set implementation
    typedef typename GridPartTraits::IndexSetType IndexSetType;

    //! \brief type of IntersectionIterator
    typedef typename GridPartTraits::IntersectionIteratorType IntersectionIteratorType;

    //! \brief type of Entity with codim=0
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = GridPartTraits :: conforming };
  public:
    //! \brief Returns const reference to the underlying grid
    const GridType & grid () const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
      return asImp().grid();
    }
    //! \brief Returns reference to the underlying grid
    GridType & grid ()
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
      return asImp().grid();
    }

    //! \brief Returns reference to index set of the underlying grid
    const IndexSetType& indexSet() const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().indexSet()));
      return asImp().indexSet();
    }

    /** \brief  Returns first iterator of the subset of the entities of codimension cd
        specified by this class
     */
    template <int cd>
    typename GridPartTraits::template Codim<cd>::IteratorType
    begin() const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template begin<cd>()));
      return asImp().template begin<cd>();
    }

    /** \brief Returns end iterator of the subset of the entities of codimension cd
        specified by this class
     */
    template <int cd>
    typename GridPartTraits::template Codim<cd>::IteratorType
    end() const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template end<cd>()));
      return asImp().template end<cd>();
    }

    //! \brief Level of the grid part
    int level() const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().level()));
      return asImp().level();
    }

    //! \brief ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().ibegin(en)));
      return asImp().ibegin(en);
    }

    //! \brief iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().iend(en)));
      return asImp().iend(en);
    }

    //! \brief corresponding communication method for grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                     InterfaceType iftype, CommunicationDirection dir) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().communicate(data,iftype,dir)));
    }

  protected:
    //! do not create explicit instances of this class
    GridPartInterface () {}

  private:
    // Barton-Nackman
    GridPartType& asImp() {
      return static_cast<GridPartType&>(*this);
    }

    // const Barton-Nackman
    const GridPartType& asImp() const {
      return static_cast<const GridPartType&>(*this);
    }
  };

  //! \brief Default implementation for the GridPart classes
  template <class GridPartTraits>
  class GridPartDefault :
    public GridPartInterface<GridPartTraits> {
  public:
    //! Grid implementation
    typedef typename GridPartTraits::GridType GridType;
    //! Index set implementation
    typedef typename GridPartTraits::IndexSetType IndexSetType;

  protected:
    //! Constructor
    GridPartDefault(GridType& grid, const IndexSetType& iset) :
      GridPartInterface<GridPartTraits>(),
      grid_(grid),
      iset_(iset) {}

    ~GridPartDefault() {}

  public:
    //! Returns const reference to the underlying grid
    const GridType& grid() const { return grid_; }

    //! Returns reference to the underlying grid
    GridType& grid() { return grid_; }

    //! Returns reference to index set of the underlying grid
    const IndexSetType& indexSet() const { return iset_; }

  private:
    GridType& grid_;
    const IndexSetType& iset_;
  };

  //! \brief Selects a specific level of a grid
  template <class GridImp, PartitionIteratorType pitype = Interior_Partition>
  class LevelGridPart :
    public GridPartDefault<LevelGridPartTraits<GridImp,pitype> > {
  public:
    //- Public typedefs and enums
    //! Corresponding type definitions
    typedef LevelGridPartTraits<GridImp,pitype> Traits;
    //! Grid implementation
    typedef typename Traits::GridType GridType;
    //! Level index set that corresponds to the grid
    typedef typename Traits::IndexSetType IndexSetType;

    //! The corresponding Intersection
    typedef typename Traits::IntersectionType IntersectionType ;

    //! The corresponding IntersectionIterator
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

    //! Struct defining the iterator types for codimension cd
    template <int cd>
    struct Codim {
      /** \brief Iterator iterating over the entities of codimension <tt>cd</tt>
          in this grid part */
      typedef typename Traits::template Codim<cd>::IteratorType IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Traits :: conforming };
  private:
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

  public:
    //- Public methods
    //! Constructor
    LevelGridPart(GridType& grid, int level ) :
      GridPartDefault<Traits>(grid,isetWrapper_),
      isetWrapper_(grid,level),
      level_(level) {}

    //! Constructor, choosing maxLevel
    LevelGridPart(const GridType& grid) :
      GridPartDefault<Traits>(grid,isetWrapper_),
      isetWrapper_(grid,grid.maxLevel()),
      level_(grid.maxLevel()) {}

    //! copy constructor
    LevelGridPart(const LevelGridPart& other) :
      GridPartDefault<Traits>(other.grid_,isetWrapper_),
      isetWrapper_(other.grid_,other.level_),
      level_(other.level_) {}

    //! Returns first iterator on this level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType begin() const {
      return this->grid().template lbegin<cd,pitype>(level_);
    }

    //! Returns end iterator on this level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType end() const {
      return this->grid().template lend<cd,pitype>(level_);
    }

    //! ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const
    {
      return en.ilevelbegin();
    }

    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const
    {
      return en.ilevelend();
    }

    //! Level which this GridPart belongs to
    int level() const { return level_; }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                     InterfaceType iftype, CommunicationDirection dir) const
    {
      this->grid().communicate(data,iftype,dir,level());
    }

  private:
    //! GridDefaultIndexSet Wrapper
    IndexSetType isetWrapper_;
    const int level_;
  };

  //! Type definitions for the LevelGridPart class
  template <class GridImp, PartitionIteratorType pitype>
  struct LevelGridPartTraits {

    /** \brief The type of the grid */
    typedef GridImp GridType;

    /** \brief The type of the corresponding grid part class */
    typedef LevelGridPart<GridImp,pitype> GridPartType;

    /** \brief The appropriate index set */
    typedef WrappedLevelIndexSet<GridType> IndexSetType;

    /** \brief The appropriate intersection */
    typedef typename GridType::Traits::
    LevelIntersection IntersectionType;

    /** \brief The appropriate intersection iterator */
    typedef typename GridType::template Codim<0>::Entity::
    LevelIntersectionIterator IntersectionIteratorType;

    /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
    template <int cd>
    struct Codim {
      /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LevelIterator IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Capabilities::isLevelwiseConforming<GridType>::v };
  };

  //! \brief Selects the leaf level of a grid
  template <class GridImp, PartitionIteratorType pitype = Interior_Partition>
  class LeafGridPart :
    public GridPartDefault<LeafGridPartTraits<GridImp,pitype> > {
  public:
    //- Public typedefs and enums
    //! Type definitions
    typedef LeafGridPartTraits<GridImp,pitype> Traits;
    //! Grid implementation type
    typedef typename Traits::GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits::IndexSetType IndexSetType;

    //! The corresponding Intersection
    typedef typename Traits::IntersectionType IntersectionType ;

    //! The corresponding IntersectionIterator
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

    //! Struct providing types of the leaf iterators on codimension cd
    template <int cd>
    struct Codim {
      //! Provide types of the leaf iterators on codimension cd
      typedef typename Traits::template Codim<cd>::IteratorType IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Traits :: conforming };
  private:
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

  public:
    //- Public methods
    //! Constructor
    LeafGridPart(GridType& grid) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

    //! copy constructor
    LeafGridPart(const LeafGridPart& other) :
      GridPartDefault<Traits>(other.grid_,isetWrapper_),
      isetWrapper_(other.grid_)
    {}

    //! Begin iterator on the leaf level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType begin() const {
      return this->grid().template leafbegin<cd,pitype>();
    }

    //! End iterator on the leaf level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType end() const {
      return this->grid().template leafend<cd,pitype>();
    }

    //! ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const
    {
      return en.ileafbegin();
    }

    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const
    {
      return en.ileafend();
    }

    //! Returns maxlevel of the grid
    int level() const { return this->grid().maxLevel(); }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                     InterfaceType iftype, CommunicationDirection dir) const
    {
      this->grid().communicate(data,iftype,dir);
    }

  private:
    //! GridDefaultIndexSet Wrapper
    IndexSetType isetWrapper_;
  };

  //! Type definitions for the LeafGridPart class
  template <class GridImp,PartitionIteratorType pitype>
  struct LeafGridPartTraits {

    /** \brief The type of the grid */
    typedef GridImp GridType;

    /** \brief The type of the corresponding grid part class */
    typedef LeafGridPart<GridImp,pitype> GridPartType;

    /** \brief The appropriate index set */
    typedef WrappedLeafIndexSet<GridType> IndexSetType;

    /** \brief The appropriate intersection */
    typedef typename GridType::Traits::
    LeafIntersection IntersectionType;

    /** \brief The appropriate intersection iterator */
    typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;

    /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
    template <int cd>
    struct Codim {
      /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
  };


  //**************************************************************
  /** \brief Selects the leaf level of a grid together with the
      HierarchicIndexSet available for ALUGrid and AlbertaGrid.
      The HierarchicIndexSet is basically the LocalIdSet of the grid
      extended by a size method to implement the IndexSet interface.
      For all other grids the default LeafIndexSet is selected.

      \deprected
   */
  template <class GridImp, PartitionIteratorType pitype = Interior_Partition>
  class HierarchicGridPart :
    public GridPartDefault< HierarchicGridPartTraits<GridImp,pitype> > {
  public:
    //- Public typedefs and enums
    //! Type definitions
    typedef HierarchicGridPartTraits<GridImp,pitype> Traits;
    //! Grid implementation type
    typedef typename Traits::GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits::IndexSetType IndexSetType;

    //! The corresponding IntersectionIterator
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

    //! Struct providing types of the leaf iterators on codimension cd
    template <int cd>
    struct Codim {
      typedef typename Traits::template Codim<cd>::IteratorType IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Traits :: conforming };
  private:
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

  public:
    //- Public methods
    //! Constructor
    HierarchicGridPart(GridType& grid) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

    //! Constructor
    HierarchicGridPart(GridType& grid, const IndexSetType & ) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

    //! copy constructor
    HierarchicGridPart(const HierarchicGridPart& other) :
      GridPartDefault<Traits>(other.grid_, isetWrapper_),
      isetWrapper_(other.grid_) {}

    //! Begin iterator on the leaf level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType begin() const {
      return this->grid().template leafbegin<cd,pitype>();
    }

    //! End iterator on the leaf level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType end() const {
      return this->grid().template leafend<cd,pitype>();
    }

    //! ibegin of corresponding intersection iterator for given entity
    IntersectionIteratorType ibegin(const EntityCodim0Type & en) const
    {
      return en.ileafbegin();
    }

    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const
    {
      return en.ileafend();
    }

    //! Returns maxlevel of the grid
    int level() const { return this->grid().maxLevel(); }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                     InterfaceType iftype, CommunicationDirection dir) const
    {
      this->grid().communicate(data,iftype,dir);
    }

  private:
    //! GridDefaultIndexSet Wrapper
    IndexSetType isetWrapper_;
  };

  //! Type definitions for the HierarchicGridPart class
  /** \deprected */
  template <class GridImp,PartitionIteratorType pitype>
  struct HierarchicGridPartTraits {
    /** \brief The type of the grid */
    typedef GridImp GridType;
    /** \brief The type of the corresponding grid part class */
    typedef HierarchicGridPart<GridImp,pitype> GridPartType;
    /** \brief The appropriate index set */
    typedef WrappedHierarchicIndexSet<GridType> IndexSetType;
    /** \brief The appropriate intersection iterator */
    typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;

    /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
    template <int cd>
    struct Codim {
      /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
  };

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
  /** @} */

} // end namespace Dune

#endif


// Error Message for dune-fem
//
// This message can only be triggered by users of dune-fem and tells them which
// file to include in the future
#ifdef DUNE_FEM_GRIDPART_HH
#error "Including <dune/grid/common/gridpart.hh> is deprecated."
#error "Include <dune/fem/gridpart/gridpart.hh> instead."
#endif
