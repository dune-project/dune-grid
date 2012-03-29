// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/entitypointer.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Traits, bool fake = Traits::fake >
    class Iterator;

    template< class Grid >
    class HierarchicIterator;



    // PartitionIteratorFilter
    // -----------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct PartitionIteratorFilter;

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Interior_Partition, Grid >
    {
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

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
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Overlap_Partition, Grid >
    {
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

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
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, All_Partition, Grid >
    {
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = All_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

      static bool apply ( const ReferenceElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Ghost_Partition, Grid >
    {
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Ghost_Partition;

      typedef typename remove_const< Grid >::type::ctype ctype;
      typedef typename remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef GenericReferenceElement< ctype, dimension > ReferenceElement;

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
      typedef typename Base::EntityImpl EntityImpl;

      using Base::hostEntityIterator_;
      using Base::entityImpl;
      using Base::grid;

    public:
      Iterator ( const Grid &grid, int level, IteratorType type )
        : Base( grid, Traits::getHostEntityIterator( grid.hostGrid(), level, type ) )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
        entityImpl() = EntityImpl( grid() );
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

    protected:
      typedef typename Base::EntityImpl EntityImpl;

      using Base::hostElementIterator_;
      using Base::entityImpl;
      using Base::grid;

    public:
      Iterator ( const Grid &grid, int level, IteratorType type )
        : Base( grid, Traits::getHostElementIterator( grid.hostGrid(), level, type ), -1 ),
          hostEnd_( Traits::getHostElementIterator( grid.hostGrid(), level, Traits::end ) ),
          hostIndexSet_( &Traits::getHostIndexSet( grid.hostGrid(), level ) )
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

        int subEntity = this->subEntity();
        while( hostElementIterator_ != hostEnd_ )
        {
          const HostElement &hostElement = *hostElementIterator_;

          const GenericReferenceElement< ctype, dimension > &refElement
            = GenericReferenceElements< ctype, dimension >::general( hostElement.type() );

          ++subEntity;
          const int count = refElement.size( codimension );
          for( ; subEntity < count; ++subEntity )
          {
            if( !Filter::apply( refElement, hostElement, subEntity ) )
              continue;

            const size_t index = hostIndexSet_->subIndex( hostElement, subEntity, codimension );
            if( !visited_[ index ] )
            {
              visited_[ index ] = true;
              entityImpl() = EntityImpl( grid(), subEntity );
              return;
            }
          }
          ++hostElementIterator_;
          subEntity = -1;
        }
        entityImpl() = EntityImpl( grid(), subEntity );
      }

    private:
      HostElementIterator hostEnd_;
      const HostIndexSet *hostIndexSet_;
      std::vector< bool > visited_;
    };



    // LeafIteratorTraits
    // ------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LeafIteratorTraits
      : public EntityPointerTraits< codim, Grid >
    {
      typedef typename EntityPointerTraits< codim, Grid >::HostGrid HostGrid;

      typedef PartitionIteratorFilter< codim, pitype, HostGrid > Filter;

      static const PartitionIteratorType Entity_Partition = pitype;
      static const PartitionIteratorType Element_Partition = Filter::Element_Partition;

      typedef typename HostGrid::template Codim< codim >
      ::template Partition< Entity_Partition >::LeafIterator
      HostEntityIterator;
      typedef typename HostGrid::template Codim< 0 >
      ::template Partition< Element_Partition >::LeafIterator
      HostElementIterator;

      typedef typename HostGrid::LeafIndexSet HostIndexSet;

      enum IteratorType { begin, end };

      static HostEntityIterator
      getHostEntityIterator ( const HostGrid &hostGrid, int level, IteratorType type )
      {
        if( type == begin )
          return hostGrid.template leafbegin< codim, Entity_Partition >();
        else
          return hostGrid.template leafend< codim, Entity_Partition >();
      }

      static HostElementIterator
      getHostElementIterator ( const HostGrid &hostGrid, int level, IteratorType type )
      {
        if( type == begin )
          return hostGrid.template leafbegin< 0, Element_Partition >();
        else
          return hostGrid.template leafend< 0, Element_Partition >();
      }

      static const HostIndexSet &getHostIndexSet ( const HostGrid &hostGrid, int level )
      {
        return hostGrid.leafIndexSet();
      }
    };



    // LevelIteratorTraits
    // -------------------

    template< int codim, PartitionIteratorType pitype, class Grid >
    struct LevelIteratorTraits
      : public EntityPointerTraits< codim, Grid >
    {
      typedef typename EntityPointerTraits< codim, Grid >::HostGrid HostGrid;

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
      getHostEntityIterator ( const HostGrid &hostGrid, int level, IteratorType type )
      {
        if( type == begin )
          return hostGrid.template lbegin< codim, Entity_Partition >( level );
        else
          return hostGrid.template lend< codim, Entity_Partition >( level );
      }

      static HostElementIterator
      getHostElementIterator ( const HostGrid &hostGrid, int level, IteratorType type )
      {
        if( type == begin )
          return hostGrid.template lbegin< 0, Element_Partition >( level );
        else
          return hostGrid.template lend< 0, Element_Partition >( level );
      }

      static const HostIndexSet &getHostIndexSet ( const HostGrid &hostGrid, int level )
      {
        return hostGrid.levelIndexSet( level );
      }
    };



    // HierarchicIteratorTraits
    // ------------------------

    template< class Grid >
    struct HierarchicIteratorTraits
      : public EntityPointerTraits< 0, Grid >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid::Traits::HierarchicIterator HostEntityIterator;
      typedef typename Traits::HostGrid::Traits::HierarchicIterator HostElementIterator;
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
      typedef typename Base::EntityImpl EntityImpl;
      typedef typename Traits::HostEntityIterator HostEntityIterator;

      using Base::hostEntityIterator_;
      using Base::entityImpl;
      using Base::grid;

    public:
      HierarchicIterator ( const Grid &grid,
                           const HostEntityIterator &hostEntityIterator )
        : Base( grid, hostEntityIterator )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
        entityImpl() = EntityImpl( grid() );
      }
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ITERATOR_HH
