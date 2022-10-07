// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_GRIDVIEW_HH
#define DUNE_ALBERTAGRID_GRIDVIEW_HH

#if HAVE_ALBERTA

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/grid/albertagrid/intersectioniterator.hh>

namespace Dune
{

  template< class GridImp >
  class AlbertaLevelGridView;

  template< class GridImp >
  class AlbertaLeafGridView;


  template< class GridImp >
  struct AlbertaLevelGridViewTraits
  {
    typedef AlbertaLevelGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid::Traits::LevelIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid::Traits::LevelIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid::Traits::LevelIntersectionIterator
    IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid::Traits::Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< cd >::template Partition< All_Partition >::LevelIterator Iterator;

      typedef typename Grid::Traits::template Codim< cd >::Entity Entity;

      typedef typename Grid::template Codim< cd >::Geometry Geometry;
      typedef typename Grid::template Codim< cd >::LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template< PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid::template Codim< cd >::template Partition< pit >::LevelIterator
        Iterator;
      };
    };

    constexpr static bool conforming = Capabilities::isLevelwiseConforming< Grid >::v;
  };


  template< class GridImp >
  class AlbertaLevelGridView
  {
    typedef AlbertaLevelGridView< GridImp > ThisType;

  public:
    typedef AlbertaLevelGridViewTraits<GridImp> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits::IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits::Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits::IntersectionIterator IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Traits::Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits::template Codim<cd> {};

    constexpr static bool conforming = Traits::conforming;

  private:
    typedef Alberta::ElementInfo< Grid::dimension > ElementInfo;

    typedef Dune::AlbertaGridLeafIntersectionIterator< GridImp > IntersectionIteratorImpl;

  public:
    AlbertaLevelGridView ( const Grid &grid, int level )
      : grid_( &grid ),
        indexSet_( &(grid.levelIndexSet( level )) ),
        level_( level )
    {}

    // use default implementation of copy constructor and assignment operator
#if 0
    AlbertaLevelGridView ( const ThisType &other )
      : grid_( other.grid_ ),
        indexSet_( other.indexSet_ ),
        level_( other.level_ )
    {}

    /** \brief assignment from other GridView on the same grid */
    ThisType &operator= ( const ThisType & other)
    {
      grid_ = other.grid_;
      indexSet_ = other.indexSet_;
      level_ = other.level_;
    }
#endif

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return *indexSet_;
    }

    /** \brief return true if current state of grid view represents a conforming grid */
    bool isConforming() const { return bool(conforming); }

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
    typename Codim< cd >::Iterator begin () const
    {
      return grid().template lbegin< cd, All_Partition >( level_ );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd >::template Partition< pit >::Iterator begin () const
    {
      return grid().template lbegin< cd, pit >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd >::Iterator end () const
    {
      return grid().template lend< cd, All_Partition >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd >::template Partition< pit >::Iterator end () const
    {
      return grid().template lend< cd, pit >( level_ );
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      if( grid().maxLevel() == 0)
      {
        typename IntersectionIteratorImpl::Begin begin;
        return IntersectionIteratorImpl( entity.impl(), begin );
      }
      else
      {
        DUNE_THROW( NotImplemented, "method ibegin not implemented on LevelGridView for AlbertaGrid." );
        typename IntersectionIteratorImpl::End end;
        return IntersectionIteratorImpl( entity.impl(), end );
      }
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      typename IntersectionIteratorImpl::End end;
      return IntersectionIteratorImpl( entity.impl(), end );
    }

    /** \brief obtain communication object */
    const Communication &comm () const
    {
      return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize ( [[maybe_unused]] int codim ) const { return 0; }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize ( [[maybe_unused]] int codim ) const { return 0; }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( [[maybe_unused]] CommDataHandleIF< DataHandleImp, DataType > &data,
                       [[maybe_unused]] InterfaceType iftype,
                       [[maybe_unused]] CommunicationDirection dir ) const
    {}

  private:
    const Grid *grid_;
    const IndexSet *indexSet_;
    int level_;
  };


  template< class GridImp >
  struct AlbertaLeafGridViewTraits
  {
    typedef AlbertaLeafGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid::Traits::LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid::Traits::LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid::Traits::LeafIntersectionIterator
    IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Grid::Traits::Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< cd >::template Partition< All_Partition >::LeafIterator
      Iterator;

      typedef typename Grid::Traits::template Codim< cd >::Entity Entity;

      typedef typename Grid::template Codim< cd >::Geometry Geometry;
      typedef typename Grid::template Codim< cd >::LocalGeometry
      LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid::template Codim< cd >::template Partition< pit >::LeafIterator
        Iterator;
      };
    };

    constexpr static bool conforming = Capabilities::isLeafwiseConforming< Grid >::v;
  };


  template< class GridImp >
  class AlbertaLeafGridView
  {
    typedef AlbertaLeafGridView< GridImp > ThisType;

  public:
    typedef AlbertaLeafGridViewTraits<GridImp> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits::IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits::Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits::IntersectionIterator IntersectionIterator;

    /** \brief type of the communication */
    typedef typename Traits::Communication Communication;

    /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits::template Codim<cd> {};

    constexpr static bool conforming = Traits::conforming;

  private:
    typedef Alberta::ElementInfo< Grid::dimension > ElementInfo;

    typedef Dune::AlbertaGridLeafIntersectionIterator< GridImp > IntersectionIteratorImpl;

  public:
    AlbertaLeafGridView ( const Grid &grid )
      : grid_( &grid ),
        indexSet_( &(grid.leafIndexSet()) )
    {}

    // use default implementation of copy constructor and assignment operator
#if 0
    AlbertaLeafGridView ( const ThisType &other )
      : grid_( other.grid_ ),
        indexSet_( other.indexSet_ )
    {}

    /** \brief assignment from other GridView on the same grid */
    ThisType &operator= ( const ThisType & other)
    {
      grid_ = other.grid_;
      indexSet_ = other.indexSet_;
    }
#endif

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return *indexSet_;
    }

    /** \brief return true if current state of grid view represents a conforming grid */
    bool isConforming() const { return bool(conforming); }

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
    typename Codim< cd >::Iterator begin () const
    {
      return grid().template leafbegin< cd, All_Partition >();
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd >::template Partition< pit >::Iterator begin () const
    {
      return grid().template leafbegin< cd, pit >();
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd >::Iterator end () const
    {
      return grid().template leafend< cd, All_Partition >();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd >::template Partition< pit >::Iterator end () const
    {
      return grid().template leafend< cd, pit >();
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      const ElementInfo elementInfo = entity.impl().elementInfo();
      assert( !!elementInfo );

#ifndef NDEBUG
      for( int i = 0; i <= Grid::dimension; ++i )
      {
        if( elementInfo.elInfo().opp_vertex[ i ] == 127 )
          DUNE_THROW( NotImplemented, "AlbertaGrid: Intersections on outside entities are not fully implemented, yet." );
      }
#endif // #ifndef NDEBUG

      typename IntersectionIteratorImpl::Begin begin;
      return IntersectionIteratorImpl( entity.impl(), begin );
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      assert( !!entity.impl().elementInfo() );
      typename IntersectionIteratorImpl::End end;
      return IntersectionIteratorImpl( entity.impl(), end );
    }

    /** \brief obtain communication object */
    const Communication &comm () const
    {
      return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize ( [[maybe_unused]] int codim ) const { return 0; }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize ( [[maybe_unused]] int codim ) const { return 0; }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( [[maybe_unused]] CommDataHandleIF< DataHandleImp, DataType > &data,
                       [[maybe_unused]] InterfaceType iftype,
                       [[maybe_unused]] CommunicationDirection dir ) const
    {}

  private:
    const Grid *grid_;
    const IndexSet *indexSet_;
  };

}

#endif // #if HAVE_ALBERTA
#endif // #ifndef DUNE_ALBERTAGRID_GRIDVIEW_HH
