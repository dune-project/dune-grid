// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_DEFAULTGRIDVIEW_HH
#define DUNE_GRID_COMMON_DEFAULTGRIDVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  template< class GridImp >
  class DefaultLevelGridView;

  template< class GridImp >
  class DefaultLeafGridView;


  template< class GridImp >
  struct DefaultLevelGridViewTraits
  {
    typedef DefaultLevelGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LevelIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LevelIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LevelIntersectionIterator
    IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid :: Traits :: Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
      :: template Codim< cd > :: template Partition< All_Partition > :: LevelIterator
      Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template< PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid :: template Codim< cd >
        :: template Partition< pit > :: LevelIterator
        Iterator;
      };
    };

    constexpr static bool conforming = Capabilities :: isLevelwiseConforming< Grid > :: v;
  };


  template< class GridImp >
  class DefaultLevelGridView
  {
    typedef DefaultLevelGridView< GridImp > ThisType;

  public:
    typedef DefaultLevelGridViewTraits<GridImp> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid :: Traits :: Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};

    constexpr static bool conforming = Traits :: conforming;

    DefaultLevelGridView ( const Grid &grid, int level )
      : grid_( &grid ),
        level_( level )
    {}

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return grid().levelIndexSet( level_ );
    }

    /** \brief return true if current state of grid view represents a conforming grid */
    bool isConforming() const { return bool(Traits::conforming); }

    /** \brief obtain number of entities in a given codimension */
    int size ( int codim ) const
    {
      return grid().size( level_, codim );
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const GeometryType &type ) const
    {
      return grid().size( level_, type );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return grid().template lbegin< cd, All_Partition >( level_ );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
    {
      return grid().template lbegin< cd, pit >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return grid().template lend< cd, All_Partition >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator end () const
    {
      return grid().template lend< cd, pit >( level_ );
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ilevelbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ilevelend();
    }

    /** \brief obtain communication object */
    const Communication &comm () const
    {
      return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
      return grid().overlapSize(level_, codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
      return grid().ghostSize(level_, codim);
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      return grid().communicate( data, iftype, dir, level_ );
    }

  private:
    const Grid *grid_;
    int level_;
  };


  template< class GridImp >
  struct DefaultLeafGridViewTraits {
    typedef DefaultLeafGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LeafIntersectionIterator
    IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid :: Traits :: Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
      :: template Codim< cd > :: template Partition< All_Partition > :: LeafIterator
      Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid :: template Codim< cd >
        :: template Partition< pit > :: LeafIterator
        Iterator;
      };
    };

    constexpr static bool conforming = Capabilities :: isLeafwiseConforming< Grid > :: v;
  };


  template< class GridImp >
  class DefaultLeafGridView
  {
    typedef DefaultLeafGridView< GridImp > ThisType;

  public:
    typedef DefaultLeafGridViewTraits<GridImp> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid :: Traits :: Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};

    constexpr static bool conforming = Traits :: conforming;

  public:
    DefaultLeafGridView ( const Grid &grid )
      : grid_( &grid )
    {}

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return grid().leafIndexSet();
    }

    /** \brief return true if current state of grid view represents a conforming grid */
    bool isConforming() const { return bool(Traits::conforming); }

    /** \brief obtain number of entities in a given codimension */
    int size ( int codim ) const
    {
      return grid().size( codim );
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const GeometryType &type ) const
    {
      return grid().size( type );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return grid().template leafbegin< cd, All_Partition >();
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
    {
      return grid().template leafbegin< cd, pit >();
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return grid().template leafend< cd, All_Partition >();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator end () const
    {
      return grid().template leafend< cd, pit >();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ileafbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ileafend();
    }

    /** \brief obtain communication object */
    const Communication &comm () const
    {
      return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
      return grid().overlapSize(codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
      return grid().ghostSize(codim);
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      return grid().communicate( data, iftype, dir );
    }

  private:
    const Grid *grid_;
  };

}

#endif
