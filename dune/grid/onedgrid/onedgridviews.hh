// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ONEDGRID_ONEDGRIDVIEWS_HH
#define DUNE_GRID_ONEDGRID_ONEDGRIDVIEWS_HH

#include <type_traits>
#include <tuple>

//- includes from dune-geometry
#include <dune/geometry/type.hh>

//- includes from dune-grid
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/grid/onedgrid/onedgridentity.hh>
#include <dune/grid/onedgrid/onedgridleafiterator.hh>

namespace Dune
{

  template< class GridImp >
  class OneDGridLevelGridView;

  template< class GridImp >
  class OneDGridLeafGridView;


  /** \brief Collect several types associated to OneDGrid LevelGridViews */
  template< class GridImp>
  struct OneDGridLevelGridViewTraits
  {
    typedef OneDGridLevelGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

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

    enum { conforming = true };
  };

  /** \brief Implementation class of LevelGridViews for OneDGrid */
  template< class GridImp >
  class OneDGridLevelGridView
  {
  public:
    typedef OneDGridLevelGridViewTraits<GridImp> Traits;

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

    OneDGridLevelGridView ( const Grid &grid, int level )
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

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return OneDGridLevelIterator<cd,All_Partition,GridImp>(const_cast<OneDEntityImp<1-cd>*>(std::get<1-cd>(grid_->entityImps_[level_]).begin()));
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
    {
      return OneDGridLevelIterator<cd,pit,GridImp>(const_cast<OneDEntityImp<1-cd>*>(std::get<1-cd>(grid_->entityImps_[level_]).begin()));
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return OneDGridLevelIterator<cd,All_Partition,GridImp>(nullptr);
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator end () const
    {
      return OneDGridLevelIterator<cd,pit,GridImp>(static_cast<OneDEntityImp<1-cd>*>(nullptr));
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

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
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
      return 0;
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {}

  private:
    const Grid *grid_;
    int level_;
  };


  /** \brief Collect several types associated to OneDGrid LeafGridViews */
  template< class GridImp>
  struct OneDGridLeafGridViewTraits
  {
    typedef OneDGridLeafGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LeafIntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
      :: template Codim< cd > :: template Partition< All_Partition > :: LeafIterator
      Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry LocalGeometry;

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

    enum { conforming = true };
  };

  /** \brief Implementation class of LeafGridViews for UGGrid */
  template< class GridImp >
  class OneDGridLeafGridView
  {
  public:
    typedef OneDGridLeafGridViewTraits<GridImp> Traits;

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

  public:
    OneDGridLeafGridView ( const Grid &grid )
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

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return OneDGridLeafIterator<cd,All_Partition,GridImp>(*grid_);
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
    {
      return OneDGridLeafIterator<cd,pit,GridImp>(*grid_);
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return OneDGridLeafIterator<cd,All_Partition,GridImp>();
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator end () const
    {
      return OneDGridLeafIterator<cd,pit,GridImp>();
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

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
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
      return 0;
    }

    /** \brief Communicate data on this view -- does nothing because OneDGrid is purely sequential */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {}

  private:
    const Grid *grid_;
  };

}

#endif
