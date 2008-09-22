// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEFAULTGRIDVIEW_HH
#define DUNE_DEFAULTGRIDVIEW_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune {

  template< class GridImp, PartitionIteratorType pitype >
  class DefaultLevelGridView;
  template< class GridImp, PartitionIteratorType pitype >
  class DefaultLeafGridView;

  template< class GridImp, PartitionIteratorType pitype >
  struct DefaultLevelGridViewTraits {
    typedef DefaultLevelGridView< GridImp, pitype > GridViewImp;

    /** \brief type of the grid */
    typedef typename remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LevelIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LevelIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LevelIntersectionIterator
    IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
      :: template Codim< cd > :: template Partition< pitype > :: LevelIterator
      Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;
      typedef typename Grid :: Traits :: template Codim< cd > :: EntityPointer
      EntityPointer;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pitype_>
      struct Partition
      {
        /** \brief The iterator needed to iterate over the entities of a given codim and
            partition type of this index set */
        typedef typename Grid :: template Codim<cd> :: template Partition<pitype_> :: LevelIterator Iterator;
      };

    };

    enum { conforming = Capabilities :: isLevelwiseConforming< Grid > :: v };
  };

  template< class GridImp, PartitionIteratorType pitype >
  class DefaultLevelGridView
  {
    typedef DefaultLevelGridView< GridImp, pitype > ThisType;

  public:
    typedef DefaultLevelGridViewTraits<GridImp,pitype> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};

    enum { conforming = Traits :: conforming };

  private:
    const Grid &grid_;
    const IndexSet &indexSet_;
    const int level_;

  public:
    DefaultLevelGridView ( const Grid &grid, int level )
      : grid_( grid ),
        indexSet_( grid.levelIndexSet( level ) ),
        level_( level )
    {}

    DefaultLevelGridView ( const ThisType &other )
      : grid_( other.grid_ ),
        indexSet_( other.indexSet_ ),
        level_( other.level_ )
    {}

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return indexSet_;
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pitype_=pitype >
    typename Codim< cd > :: template Partition< pitype_ > :: Iterator begin () const
    {
      return grid().template lbegin< cd, pitype_ >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pitype_=pitype >
    typename Codim< cd > :: template Partition< pitype_ > :: Iterator end () const
    {
      return grid().template lend< cd, pitype_ >( level_ );
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.ilevelbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.ilevelend();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
    {
      return grid().comm();
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      return grid().communicate( data, iftype, dir, level_ );
    }
  };


  template< class GridImp, PartitionIteratorType pitype >
  struct DefaultLeafGridViewTraits {
    typedef DefaultLeafGridView< GridImp, pitype > GridViewImp;

    /** \brief type of the grid */
    typedef typename remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LeafIntersectionIterator
    IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
      :: template Codim< cd > :: template Partition< pitype > :: LeafIterator
      Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;
      typedef typename Grid :: Traits :: template Codim< cd > :: EntityPointer
      EntityPointer;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pitype_>
      struct Partition
      {
        /** \brief The iterator needed to iterate over the entities of a given codim and
            partition type of this index set */
        typedef typename Grid :: template Codim<cd> :: template Partition<pitype_> :: LeafIterator Iterator;
      };

    };

    enum { conforming = Capabilities :: isLeafwiseConforming< Grid > :: v };
  };

  template< class GridImp, PartitionIteratorType pitype >
  class DefaultLeafGridView
  {
    typedef DefaultLeafGridView< GridImp, pitype > ThisType;

  public:
    typedef DefaultLeafGridViewTraits<GridImp,pitype> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};

    enum { conforming = Traits :: conforming };

  private:
    const Grid &grid_;
    const IndexSet &indexSet_;

  public:
    DefaultLeafGridView ( const Grid &grid )
      : grid_( grid ),
        indexSet_( grid.leafIndexSet() )
    {}

    DefaultLeafGridView ( const ThisType &other )
      : grid_( other.grid_ ),
        indexSet_( other.indexSet_ )
    {}

  private:
    // prohibit assignment
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return indexSet_;
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pitype_=pitype >
    typename Codim< cd > :: template Partition< pitype_ > :: Iterator begin () const
    {
      return grid().template leafbegin< cd, pitype_ >();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pitype_=pitype >
    typename Codim< cd > :: template Partition< pitype_ > :: Iterator end () const
    {
      return grid().template leafend< cd, pitype_ >();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.ileafbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.ileafend();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
    {
      return grid().comm();
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      return grid().communicate( data, iftype, dir );
    }
  };

}

#endif
