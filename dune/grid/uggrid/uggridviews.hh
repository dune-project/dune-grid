// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UGGRID_UGGRIDVIEWS_HH
#define DUNE_GRID_UGGRID_UGGRIDVIEWS_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

    template< class GridImp>
    class UGGridLevelGridView;

    template< class GridImp>
    class UGGridLeafGridView;


    /** \brief Collect several types associated to UGGrid LevelGridViews */
    template< class GridImp>
    struct UGGridLevelGridViewTraits
    {
      typedef UGGridLevelGridView< GridImp > GridViewImp;

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


    /** \brief Implementation class of LevelGridViews for UGGrid */
    template< class GridImp>
    class UGGridLevelGridView
    {
      typedef UGGridLevelGridView<GridImp> This;

    public:
      typedef UGGridLevelGridViewTraits<GridImp> Traits;

      /** \brief type of the grid */
      typedef typename Traits::Grid Grid;

      /** \brief type of the index set */
      typedef typename Traits :: IndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Traits :: Intersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Traits :: IntersectionIterator IntersectionIterator;

      /** \brief type of the communication */
      typedef typename Traits :: Communication Communication;

      /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
      [[deprecated("Use Communication instead!")]]
      typedef Communication CollectiveCommunication;

      /** \brief Codim Structure */
      template< int cd >
      struct Codim : public Traits :: template Codim<cd> {};

      constexpr static bool conforming = Traits :: conforming;

      static constexpr int dimension = Grid::dimension;

      UGGridLevelGridView ( const Grid &grid, int level )
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

      bool isConforming() const
      {
        return Traits::conforming;
      }

      /** \brief obtain begin iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator begin () const
      {
        if (!grid().multigrid_)
          DUNE_THROW(GridError, "The grid has not been properly initialized!");

        if (!grid().multigrid_->grids[level_])
          DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level_ << " requested!");

        return UGGridLevelIterator<cd, All_Partition, const Grid>( grid(), level_ );
      }

      /** \brief obtain begin iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
      {
        if (!grid().multigrid_)
          DUNE_THROW(GridError, "The grid has not been properly initialized!");

        if (!grid().multigrid_->grids[level_])
          DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level_ << " requested!");

        return UGGridLevelIterator<cd, pit, const Grid>( grid(), level_ );
      }

      /** \brief obtain end iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator end () const
      {
        return UGGridLevelIterator<cd, All_Partition, const Grid>();
      }

      /** \brief obtain end iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator end () const
      {
        return UGGridLevelIterator<cd, pit, const Grid>();
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
        return 0;
      }

      /** \brief Return size of the ghost region for a given codim on the grid view.  */
      int ghostSize(int codim) const
      {
        return (codim==0) ? 1 : 0;
      }

      /** communicate data on this view */
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &dataHandle, InterfaceType iftype, CommunicationDirection dir ) const
      {
#ifdef ModelP
        Hybrid::forEach(std::make_index_sequence< Grid::dimension+1 >{},
          [&](auto codim)
        {
          if (dataHandle.contains(Grid::dimension, codim))
            grid().template communicateUG_<codim>(*this, level_, dataHandle, iftype, dir);
        });
#endif // ModelP
      }

    private:
      const Grid *grid_;
      int level_;
    };


    template< class GridImp>
    struct UGGridLeafGridViewTraits {
      typedef UGGridLeafGridView< GridImp > GridViewImp;

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


    /** \brief Implementation class of LeafGridViews for UGGrid */
    template< class GridImp >
    class UGGridLeafGridView
    {
      typedef UGGridLeafGridView<GridImp> This;

    public:
      typedef UGGridLeafGridViewTraits<GridImp> Traits;

      /** \brief type of the grid */
      typedef typename Traits::Grid Grid;

      /** \brief type of the index set */
      typedef typename Traits :: IndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Traits :: Intersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Traits :: IntersectionIterator IntersectionIterator;

      /** \brief type of the communication */
      typedef typename Traits :: Communication Communication;

      /** \deprecated Use Communication instead! Will be removed after Dune 2.9. */
      [[deprecated("Use Communication instead!")]]
      typedef Communication CollectiveCommunication;

      /** \brief Codim Structure */
      template< int cd >
      struct Codim : public Traits :: template Codim<cd> {};

      constexpr static bool conforming = Traits :: conforming;

      static constexpr int dimension = Grid::dimension;

    public:
      UGGridLeafGridView ( const Grid &grid )
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

      bool isConforming() const
      {
        return Traits::conforming;
      }

      /** \brief obtain begin iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator begin () const
      {
        return UGGridLeafIterator<cd, All_Partition, const Grid>(grid());
      }

      /** \brief obtain begin iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
      {
        return UGGridLeafIterator<cd, pit, const Grid>(grid());
      }

      /** \brief obtain end iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator end () const
      {
        return UGGridLeafIterator<cd, All_Partition, const Grid>();
      }

      /** \brief obtain end iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator end () const
      {
        return UGGridLeafIterator<cd, pit, const Grid>();
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
        return 0;
      }

      /** \brief Return size of the ghost region for a given codim on the grid view.  */
      int ghostSize(int codim) const
      {
        return (codim==0) ? 1 : 0;
      }

      /** communicate data on this view */
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &dataHandle, InterfaceType iftype, CommunicationDirection dir ) const
      {
#ifdef ModelP
        Hybrid::forEach(std::make_index_sequence< Grid::dimension+1 >{},
          [&](auto codim)
        {
          if (dataHandle.contains(Grid::dimension, codim))
            grid().template communicateUG_<codim>(*this, -1, dataHandle, iftype, dir);
        });
#endif // ModelP
      }

    private:
      const Grid *grid_;
    };

}

#endif // #ifndef DUNE_GRID_UGGRID_UGGRIDVIEWS_HH
