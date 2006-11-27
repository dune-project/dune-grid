// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDPART_HH
#define DUNE_GRIDPART_HH

//- System includes

//- Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/common/bartonnackmanifcheck.hh>

/** @file
   @author Robert Kloefkorn
   @brief Provides views of grid via grid parts, heavily used in the
   dune-fem module.
 */

namespace Dune {

  // Forward declarations
  template <class GridImp, PartitionIteratorType pitype>
  class LevelGridPartTraits;
  template <class GridImp, PartitionIteratorType pitype>
  class LeafGridPartTraits;
  template <class GridImp, PartitionIteratorType pitype>
  class HierarchicGridPartTraits;
  template <class GridImp,class IndexSetImp, PartitionIteratorType pitype>
  struct DefaultGridPartTraits;

  //! \brief Interface for the GridPart classes
  //! A GridPart class allows to access only a specific subset of a grid's
  //! entities. A GridPart implementation provides the corresponding index set
  //! and a begin/end iterator pair for accessing those entities, the
  //! corresponding intersection iterators and a appropriate communication
  //! method.
  //! GridParts are used to parametrize spaces (see DiscreteFunctionSpaceDefault).
  template <class GridPartTraits>
  class GridPartInterface {
  public:
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
    //! \brief Returns reference to the underlying grid
    const GridType & grid () const
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
      CHECK_INTERFACE_IMPLEMENTATION((asImp().begin()));
      return asImp().begin();
    }

    /** \brief Returns end iterator of the subset of the entities of codimension cd
        specified by this class
     */
    template <int cd>
    typename GridPartTraits::template Codim<cd>::IteratorType
    end() const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().end()));
      return asImp().end();
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
    //! do not create explict instances of this class
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
    GridPartDefault(const GridType& grid, const IndexSetType& iset) :
      GridPartInterface<GridPartTraits>(),
      grid_(grid),
      iset_(iset) {}

    ~GridPartDefault() {}

  public:
    //! Returns reference to the underlying grid
    const GridType& grid() const { return grid_; }

    //! Returns reference to index set of the underlying grid
    const IndexSetType& indexSet() const { return iset_; }

  private:
    const GridType& grid_;
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

    //! The corresponding IntersectionIterator
    typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

    //! Struct defining the iterator types for codimension cd
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
    LevelGridPart(const GridType& grid, int level ) :
      GridPartDefault<Traits>(grid,isetWrapper_),
      isetWrapper_(grid,level),
      level_(level) {}

    //! Constructor, chosing maxLevel
    LevelGridPart(const GridType& grid) :
      GridPartDefault<Traits>(grid,isetWrapper_),
      isetWrapper_(grid,grid.maxLevel()),
      level_(grid.maxLevel()) {}

    //! Returns first iterator on a given level
    template <int cd>
    typename Traits::template Codim<cd>::IteratorType begin() const {
      return this->grid().template lbegin<cd,pitype>(level_);
    }

    //! Returns end iterator on a given level
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
    typedef GridImp GridType;
    typedef LevelGridPart<GridImp,pitype> GridPartType;
    typedef WrappedLevelIndexSet<GridType> IndexSetType;
    typedef typename GridType::template Codim<0>::Entity::
    LevelIntersectionIterator IntersectionIteratorType;

    template <int cd>
    struct Codim {
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
    LeafGridPart(const GridType& grid) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

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
    typedef GridImp GridType;
    typedef LeafGridPart<GridImp,pitype> GridPartType;
    typedef WrappedLeafIndexSet<GridType> IndexSetType;
    typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;

    template <int cd>
    struct Codim {
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
  };


  //**************************************************************
  //! \brief Selects the leaf level of a grid
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
    HierarchicGridPart(const GridType& grid) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

    //! Constructor
    HierarchicGridPart(const GridType& grid, const IndexSetType & ) :
      GridPartDefault<Traits>(grid, isetWrapper_),
      isetWrapper_(grid) {}

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
  struct HierarchicGridPartTraits {
    typedef GridImp GridType;
    typedef HierarchicGridPart<GridImp,pitype> GridPartType;
    typedef WrappedHierarchicIndexSet<GridType> IndexSetType;
    typedef typename GridType::template Codim<0>::Entity::
    LeafIntersectionIterator IntersectionIteratorType;

    template <int cd>
    struct Codim {
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections
    enum { conforming = Capabilities::isLeafwiseConforming<GridType>::v };
  };


  //! quick hack, to be revised by me
  //! \brief Selects the leaf level of a grid
  template <class GridImp, class IndexSetImp , PartitionIteratorType pitype = Interior_Partition>
  class DefaultGridPart :
    public GridPartDefault<DefaultGridPartTraits<GridImp,IndexSetImp,pitype> > {
  public:
    //- Public typedefs and enums
    //! Type definitions
    typedef DefaultGridPartTraits<GridImp,IndexSetImp,pitype> Traits;
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
  private:
    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

  public:
    //- Public methods
    //! Constructor
    DefaultGridPart(const GridType& grid, const IndexSetType & iset ) DUNE_DEPRECATED
      : GridPartDefault<Traits>(grid, iset) {}

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
      return en.ibegin();
    }

    //! iend of corresponding intersection iterator for given entity
    IntersectionIteratorType iend(const EntityCodim0Type & en) const
    {
      return en.iend();
    }

    //! Level of the grid part
    int level() const { return this->grid().maxLevel(); }

    //! corresponding communication method for this grid part
    template <class DataHandleImp,class DataType>
    void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                     InterfaceType iftype, CommunicationDirection dir) const
    {
      this->grid().communicate(data,iftype,dir);
    }
  };

  //! Type definitions for the LeafGridPart class
  template <class GridImp,class IndexSetImp, PartitionIteratorType pitype>
  struct DefaultGridPartTraits {
    typedef GridImp GridType;
    typedef DefaultGridPart<GridImp,IndexSetImp,pitype> GridPartType;
    typedef IndexSetImp IndexSetType;
    typedef typename GridType::template Codim<0>::Entity::
    IntersectionIterator IntersectionIteratorType;

    template <int cd>
    struct Codim {
      typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator IteratorType;
    };
  };

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION

} // end namespace Dune

#endif
