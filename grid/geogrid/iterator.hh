// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/grid/geogrid/entitypointer.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Traits, bool fake = Traits::fake >
    class Iterator;

    template< int codim, PartitionIteratorType pitype, class Grid >
    class LeafIterator;

    template< int codim, PartitionIteratorType pitype, class Grid >
    class LevelIterator;

    template< class Grid >
    class HierarchicIterator;



    // PartitionIteratorFilter
    // -----------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct PartitionIteratorFilter;

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Interior_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j )->partitionType();
          if( type == InteriorEntity )
            return true;
        }
        return false;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, InteriorBorder_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Overlap_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        if( element.partitionType() == InteriorEntity )
          return true;

        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j )->partitionType();
          if( (type == OverlapEntity) || (type == BorderEntity) )
            return true;
        }
        return false;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, OverlapFront_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, All_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = All_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Ghost_Partition, Grid >
    {
      static const int dimension = Grid::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Ghost_Partition;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< typename Grid::ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j )->partitionType();
          if( type == GhostEntity )
            return true;
        }
        return false;
      }
    };



    // Iterator (real)
    // ---------------

    template< class Traits >
    class Iterator< Traits, false >
      : public EntityPointer< Traits, false >
    {
      typedef EntityPointer< Traits, false > Base;

      typedef typename Traits::Grid Grid;

    public:
      typedef typename Traits::IteratorType IteratorType;

    protected:
      using Base::hostEntityIterator_;
      using Base::releaseEntity;

    public:
      Iterator ( const Grid &grid, int level, IteratorType type )
        : Base( grid, Traits :: getHostEntityIterator( grid, level, type ) )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
        releaseEntity();
      }
    };



    // Iterator (fake)
    // ---------------

    template< class Traits >
    class Iterator< Traits, true >
      : public EntityPointer< Traits, true >
    {
      typedef EntityPointer< Traits, true > Base;

      typedef typename Traits::Grid Grid;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;

      typedef typename Traits::IteratorType IteratorType;

    private:
      typedef typename Traits::Filter Filter;

      typedef typename Traits::HostElement HostElement;
      typedef typename Traits::HostElementIterator HostElementIterator;
      typedef typename Traits::HostIndexSet HostIndexSet;

      HostElementIterator hostEnd_;
      const HostIndexSet *hostIndexSet_;
      std::vector< bool > visited_;

    protected:
      using Base::hostElementIterator_;
      using Base::subEntity_;
      using Base::releaseEntity;

    public:
      Iterator ( const Grid &grid, int level, IteratorType type )
        : Base( grid, Traits::getHostElementIterator( grid, level, type ), -1 ),
          hostEnd_( Traits::getHostElementIterator( grid, level, Traits::end ) ),
          hostIndexSet_( &Traits::getHostIndexSet( grid, level ) )
      {
        if( hostElementIterator_ != hostEnd_ )
        {
          visited_.resize( hostIndexSet_->size( codimension ), false );
          increment();
        }
      }

      void increment ()
      {
        typedef typename Traits::ctype ctype;

        while( hostElementIterator_ != hostEnd_ )
        {
          const HostElement &hostElement = *hostElementIterator_;

          const GenericReferenceElement< ctype, dimension > &refElement
            = GenericReferenceElements< ctype, dimension >::general( hostElement.type() );

          ++subEntity_;
          const int count = refElement.size( codimension );
          for( ; subEntity_ < count; ++subEntity_ )
          {
            if( !Filter::apply( refElement, hostElement, subEntity_ ) )
              continue;

            const size_t index = hostIndexSet_->subIndex( hostElement, subEntity_, codimension );
            if( !visited_[ index ] )
            {
              visited_[ index ] = true;
              releaseEntity();
              return;
            }
          }
          ++hostElementIterator_;
          subEntity_ = -1;
        }
        releaseEntity();
      }
    };



    // LeafIteratorTraits
    // ------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LeafIteratorTraits
      : public EntityPointerTraits< codim, Grid >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

      typedef PartitionIteratorFilter< codim, pitype, HostGrid > Filter;

      static const PartitionIteratorType Entity_Partition = pitype;
      static const PartitionIteratorType Element_Partition = Filter::Element_Partition;

      typedef typename HostGrid::template Codim< codim >
      ::template Partition< Entity_Partition >::LeafIterator
      HostEntityIterator;
      typedef typename HostGrid :: template Codim< 0 >
      ::template Partition< Element_Partition >::LeafIterator
      HostElementIterator;

      typedef typename HostGrid::LeafIndexSet HostIndexSet;

      enum IteratorType { begin, end };

      static HostEntityIterator
      getHostEntityIterator ( const Grid &grid, int level, IteratorType type )
      {
        if( type == begin )
          return grid.hostGrid().template leafbegin< codim, Entity_Partition >();
        else
          return grid.hostGrid().template leafend< codim, Entity_Partition >();
      }

      static HostElementIterator
      getHostElementIterator ( const Grid &grid, int level, IteratorType type )
      {
        if( type == begin )
          return grid.hostGrid().template leafbegin< 0, Element_Partition >();
        else
          return grid.hostGrid().template leafend< 0, Element_Partition >();
      }

      static const HostIndexSet &getHostIndexSet ( const Grid &grid, int level )
      {
        return grid.hostGrid().leafIndexSet();
      }
    };



    // GridLeafIterator
    // ----------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LeafIterator
      : public Iterator< LeafIteratorTraits< codim, pitype, Grid > >
    {
      typedef LeafIteratorTraits< codim, pitype, Grid > Traits;

      typedef typename Traits::IteratorType IteratorType;

      LeafIterator ( const Grid &grid, IteratorType type )
        : Iterator< Traits >( grid, -1, type )
      {}
    };



    // LevelIteratorTraits
    // -------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LevelIteratorTraits
      : public EntityPointerTraits< codim, Grid >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

      typedef PartitionIteratorFilter< codim, pitype, HostGrid > Filter;

      static const PartitionIteratorType Entity_Partition = pitype;
      static const PartitionIteratorType Element_Partition = Filter::Element_Partition;

      typedef typename HostGrid::template Codim< codim >
      ::template Partition< Entity_Partition >::LevelIterator
      HostEntityIterator;
      typedef typename HostGrid::template Codim< 0 >
      ::template Partition< Element_Partition >::LevelIterator
      HostElementIterator;

      typedef typename HostGrid::LevelIndexSet HostIndexSet;

      enum IteratorType { begin, end };

      static HostEntityIterator
      getHostEntityIterator ( const Grid &grid, int level, IteratorType type )
      {
        if( type == begin )
          return grid.hostGrid().template lbegin< codim, Entity_Partition >( level );
        else
          return grid.hostGrid().template lend< codim, Entity_Partition >( level );
      }

      static HostElementIterator
      getHostElementIterator ( const Grid &grid, int level, IteratorType type )
      {
        if( type == begin )
          return grid.hostGrid().template lbegin< 0, Element_Partition >( level );
        else
          return grid.hostGrid().template lend< 0, Element_Partition >( level );
      }

      static const HostIndexSet &getHostIndexSet ( const Grid &grid, int level )
      {
        return grid.hostGrid().levelIndexSet( level );
      }
    };



    // LevelIterator
    // -------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LevelIterator
      : public Iterator< LevelIteratorTraits< codim, pitype, Grid > >
    {
      typedef LevelIteratorTraits< codim, pitype, Grid > Traits;

      typedef typename Traits::IteratorType IteratorType;

      LevelIterator ( const Grid &grid, int level, IteratorType type )
        : Iterator< Traits >( grid, level, type )
      {}
    };



    // HierarchicIteratorTraits
    // ------------------------

    template< class Grid >
    struct HierarchicIteratorTraits
      : public EntityPointerTraits< 0, Grid >
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

      typedef typename Traits :: HostGrid :: Traits :: HierarchicIterator
      HostEntityIterator;
      typedef typename Traits :: HostGrid :: Traits :: HierarchicIterator
      HostElementIterator;
    };



    // HierarchicIterator
    // ------------------

    template< class Grid >
    class HierarchicIterator
      : public EntityPointer< HierarchicIteratorTraits< Grid > >
    {
      typedef HierarchicIteratorTraits< Grid > Traits;

      typedef EntityPointer< Traits > Base;

    protected:
      typedef typename Traits :: HostEntityIterator HostEntityIterator;

      using Base :: hostEntityIterator_;
      using Base :: releaseEntity;

    public:
      HierarchicIterator ( const Grid &grid,
                           const HostEntityIterator &hostEntityIterator )
        : Base( grid, hostEntityIterator )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
        releaseEntity();
      }
    };

  }

}

#endif
