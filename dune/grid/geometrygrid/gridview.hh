// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRIDVIEW_HH
#define DUNE_GEOGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridImp, PartitionIteratorType pitype >
    class LevelGridView;

    template< class GridImp, PartitionIteratorType pitype >
    class LeafGridView;



    // LevelGridViewTraits
    // -------------------

    template< class GridImp, PartitionIteratorType pitype >
    struct LevelGridViewTraits
    {
      typedef LevelGridView< GridImp, pitype > GridViewImp;

      typedef typename remove_const< GridImp >::type Grid;

      typedef typename Grid::Traits::LevelIndexSet IndexSet;

      typedef typename Grid::Traits::LevelIntersection Intersection;

      typedef typename Grid::Traits::LevelIntersectionIterator
      IntersectionIterator;

      typedef typename Grid::Traits::CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim
      {
        typedef typename Grid::Traits
        ::template Codim< cd >::template Partition< pitype >::LevelIterator
        Iterator;

        typedef typename Grid::Traits::template Codim< cd >::Entity Entity;
        typedef typename Grid::Traits::template Codim< cd >::EntityPointer
        EntityPointer;

        typedef typename Grid::template Codim< cd >::Geometry Geometry;
        typedef typename Grid::template Codim< cd >::LocalGeometry
        LocalGeometry;

        template< PartitionIteratorType pit >
        struct Partition
        {
          typedef typename Grid::template Codim< cd >
          ::template Partition< pit >::LevelIterator
          Iterator;
        };
      };

      static const bool conforming = Capabilities::isLevelwiseConforming< Grid >::v;
    };



    // LevelGridView
    // -------------

    template< class GridImp, PartitionIteratorType pitype >
    class LevelGridView
    {
      typedef LevelGridView< GridImp, pitype > This;

    public:
      typedef LevelGridViewTraits< GridImp, pitype > Traits;

      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;

      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim : public Traits::template Codim<cd> {};

      static const bool conforming = Traits::conforming;

      LevelGridView ( const Grid &grid, int level )
        : grid_( &grid ),
          level_( level )
      {}

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        return grid().levelIndexSet( level_ );
      }

      int size ( int codim ) const
      {
        return grid().size( level_, codim );
      }

      int size ( const GeometryType &type ) const
      {
        return grid().size( level_, type );
      }

      template< int cd >
      typename Codim< cd >::Iterator begin () const
      {
        return grid().template lbegin< cd, pitype >( level_ );
      }

      template< int cd, PartitionIteratorType pit >
      typename Codim< cd >::template Partition< pit >::Iterator begin () const
      {
        return grid().template lbegin< cd, pit >( level_ );
      }

      template< int cd >
      typename Codim< cd >::Iterator end () const
      {
        return grid().template lend< cd, pitype >( level_ );
      }

      template< int cd, PartitionIteratorType pit >
      typename Codim< cd >::template Partition< pit >::Iterator end () const
      {
        return grid().template lend< cd, pit >( level_ );
      }

      IntersectionIterator
      ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        return entity.ilevelbegin();
      }

      IntersectionIterator
      iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        return entity.ilevelend();
      }

      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }

      int overlapSize ( int codim ) const
      {
        return grid().overlapSize( level_, codim );
      }

      int ghostSize ( int codim ) const
      {
        return grid().ghostSize( level_, codim );
      }

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



    // LeafGridViewTraits
    // ------------------

    template< class GridImp, PartitionIteratorType pitype >
    struct LeafGridViewTraits
    {
      typedef LeafGridView< GridImp, pitype > GridViewImp;

      typedef typename remove_const<GridImp>::type Grid;

      typedef typename Grid::Traits::LeafIndexSet IndexSet;

      typedef typename Grid::Traits::LeafIntersection Intersection;

      typedef typename Grid::Traits::LeafIntersectionIterator
      IntersectionIterator;

      typedef typename Grid::Traits::CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim
      {
        typedef typename Grid::Traits
        ::template Codim< cd >::template Partition< pitype >::LeafIterator
        Iterator;

        typedef typename Grid::Traits::template Codim< cd >::Entity Entity;
        typedef typename Grid::Traits::template Codim< cd >::EntityPointer
        EntityPointer;

        typedef typename Grid::template Codim< cd >::Geometry Geometry;
        typedef typename Grid::template Codim< cd >::LocalGeometry
        LocalGeometry;

        template <PartitionIteratorType pit >
        struct Partition
        {
          typedef typename Grid::template Codim< cd >
          ::template Partition< pit >::LeafIterator
          Iterator;
        };
      };

      static const bool conforming = Capabilities::isLeafwiseConforming< Grid >::v;
    };



    // LeafGridView
    // ------------

    template< class GridImp, PartitionIteratorType pitype >
    class LeafGridView
    {
      typedef LeafGridView< GridImp, pitype > This;

    public:
      typedef LeafGridViewTraits<GridImp,pitype> Traits;

      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;

      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim : public Traits::template Codim<cd> {};

      static const bool conforming = Traits::conforming;

    public:
      LeafGridView ( const Grid &grid )
        : grid_( &grid )
      {}

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        return grid().leafIndexSet();
      }

      int size ( int codim ) const
      {
        return grid().size( codim );
      }

      int size ( const GeometryType &type ) const
      {
        return grid().size( type );
      }

      template< int cd >
      typename Codim< cd >::Iterator begin () const
      {
        return grid().template leafbegin< cd, pitype >();
      }

      template< int cd, PartitionIteratorType pit >
      typename Codim< cd >::template Partition< pit >::Iterator begin () const
      {
        return grid().template leafbegin< cd, pit >();
      }

      template< int cd >
      typename Codim< cd >::Iterator end () const
      {
        return grid().template leafend< cd, pitype >();
      }

      template< int cd, PartitionIteratorType pit >
      typename Codim< cd >::template Partition< pit >::Iterator end () const
      {
        return grid().template leafend< cd, pit >();
      }

      IntersectionIterator
      ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        return entity.ileafbegin();
      }

      IntersectionIterator
      iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        return entity.ileafend();
      }

      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }

      int overlapSize ( int codim ) const
      {
        return grid().overlapSize( codim );
      }

      int ghostSize ( int codim ) const
      {
        return grid().ghostSize( codim );
      }

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

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GRIDVIEW_HH
