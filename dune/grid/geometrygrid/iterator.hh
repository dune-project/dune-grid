// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

#include <cassert>

#include <type_traits>
#include <utility>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/entity.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class HostGridView, int codim, PartitionIteratorType pitype, class Grid,
              bool fake = !Capabilities::hasHostEntity< Grid, codim >::v >
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
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j ).partitionType();
          if( type == InteriorEntity )
            return true;
        }
        return false;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, InteriorBorder_Partition, Grid >
    {
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Interior_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Overlap_Partition, Grid >
    {
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        if( element.partitionType() == InteriorEntity )
          return true;

        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j ).partitionType();
          if( (type == OverlapEntity) || (type == BorderEntity) )
            return true;
        }
        return false;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, OverlapFront_Partition, Grid >
    {
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Overlap_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, All_Partition, Grid >
    {
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = All_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        return true;
      }
    };

    template< int codim, class Grid >
    struct PartitionIteratorFilter< codim, Ghost_Partition, Grid >
    {
      static const int dimension = std::remove_const< Grid >::type::dimension;
      static const int codimension = codim;

      static const PartitionIteratorType Element_Partition = Ghost_Partition;

      typedef typename std::remove_const< Grid >::type::ctype ctype;
      typedef typename std::remove_const< Grid >::type::Traits::template Codim< 0 >::Entity Element;
      typedef Transitional::ReferenceElement< ctype, Dim<dimension> > RefElement;

      static bool apply ( const RefElement &refElement,
                          const Element &element, int subEntity )
      {
        const int size = refElement.size( subEntity, codim, dimension );
        for( int i = 0; i < size; ++i )
        {
          const int j = refElement.subEntity( subEntity, codim, i, dimension );
          PartitionType type = element.template subEntity< dimension >( j ).partitionType();
          if( type == GhostEntity )
            return true;
        }
        return false;
      }
    };



    // Iterator (real)
    // ---------------

    template< class HostGridView, int codim, PartitionIteratorType pitype, class G >
    class Iterator< HostGridView, codim, pitype, G, false >
    {
      typedef typename std::remove_const< G >::type::Traits Traits;

    public:
      typedef typename Traits::Grid Grid;

      static const int codimension = codim;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

      static const bool fake = false;

    private:
      typedef GeoGrid::Entity< codimension, Traits::dimension, G > EntityImpl;

      typedef typename HostGridView::template Codim< codim >::template Partition< pitype >::Iterator HostEntityIterator;

    public:
      Iterator () : grid_( nullptr ) {}

      Iterator ( const Grid &grid, HostEntityIterator hostEntityIterator )
        : grid_( &grid ),
          hostEntityIterator_( std::move( hostEntityIterator ) )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
      }

      bool equals ( const Iterator &rhs ) const
      {
        return hostEntityIterator_ == rhs.hostEntityIterator_;
      }

      Entity dereference () const
      {
        return EntityImpl( grid(), *hostEntityIterator_ );
      }

      int level () const { return hostEntityIterator_.level(); }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      static Iterator begin ( const Grid &grid, const HostGridView &hostGridView )
      {
        HostEntityIterator hostEntityIterator = hostGridView.template begin< codimension, pitype >();
        return Iterator( grid, std::move( hostEntityIterator ) );
      }

      static Iterator end ( const Grid &grid, const HostGridView &hostGridView )
      {
        HostEntityIterator hostEntityIterator = hostGridView.template end< codimension, pitype >();
        return Iterator( grid, std::move( hostEntityIterator ) );
      }

    private:
      const Grid *grid_;
      HostEntityIterator hostEntityIterator_;
    };



    // Iterator (fake)
    // ---------------

    template< class HostGridView, int codim, PartitionIteratorType pitype, class G >
    class Iterator< HostGridView, codim, pitype, G, true >
    {
      typedef typename std::remove_const< G >::type::Traits Traits;

    public:
      typedef typename Traits::Grid Grid;

      static const int codimension = codim;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

    private:
      typedef GeoGrid::Entity< codimension, Traits::dimension, G > EntityImpl;

      typedef PartitionIteratorFilter< codim, pitype, typename HostGridView::Grid > Filter;

      typedef typename HostGridView::template Codim<0>::template Partition< Filter::Element_Partition >::Iterator HostElementIterator;
      typedef typename HostElementIterator::Entity HostElement;
      typedef typename HostGridView::IndexSet HostIndexSet;

    public:
      Iterator () : grid_( nullptr ), subEntity_( -1 ), hostIndexSet_( nullptr ) {}

      Iterator ( const Grid &grid, HostElementIterator hostElementIterator, HostElementIterator hostEnd, const HostIndexSet &hostIndexSet )
        : grid_( &grid ),
          hostElementIterator_( hostElementIterator ),
          hostEnd_( hostEnd ),
          subEntity_( -1 ),
          hostIndexSet_( &hostIndexSet )
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

          auto refElement = referenceElement< ctype, Traits::dimension >( hostElement.type() );

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
              return;
            }
          }
          ++hostElementIterator_;
          subEntity_ = -1;
        }
      }

      bool equals ( const Iterator &rhs ) const
      {
        return hostElementIterator_ == rhs.hostElementIterator_ && ( hostElementIterator_ == hostEnd_ || subEntity_ == rhs.subEntity_ );
      }

      Entity dereference () const
      {
        return EntityImpl( grid(), *hostElementIterator_, subEntity_ );
      }

      int level () const { return hostElementIterator_.level(); }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      static Iterator begin ( const Grid &grid, const HostGridView &hostGridView )
      {
        HostElementIterator first = hostGridView.template begin< 0, Filter::Element_Partition >();
        HostElementIterator last = hostGridView.template end< 0, Filter::Element_Partition >();
        return Iterator( grid, std::move( first ), std::move( last ), hostGridView.indexSet() );
      }

      static Iterator end ( const Grid &grid, const HostGridView &hostGridView )
      {
        HostElementIterator last = hostGridView.template end< 0, Filter::Element_Partition >();
        return Iterator( grid, last, last, hostGridView.indexSet() );
      }

    private:
      const Grid *grid_;
      HostElementIterator hostElementIterator_, hostEnd_;
      int subEntity_;
      const HostIndexSet *hostIndexSet_;
      std::vector< bool > visited_;
    };



    // HierarchicIterator
    // ------------------

    template< class G >
    class HierarchicIterator
    {
      typedef typename std::remove_const< G >::type::Traits Traits;

    public:
      typedef typename Traits::Grid Grid;

      static const int codimension = 0;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

    private:
      typedef GeoGrid::Entity< codimension, Traits::dimension, G > EntityImpl;

      typedef typename Grid::HostGrid::HierarchicIterator HostEntityIterator;

    public:
      HierarchicIterator () : grid_( nullptr ) {}

      HierarchicIterator ( const Grid &grid, HostEntityIterator hostEntityIterator )
        : grid_( &grid ),
          hostEntityIterator_( std::move( hostEntityIterator ) )
      {}

      void increment ()
      {
        ++hostEntityIterator_;
      }

      bool equals ( const HierarchicIterator &rhs ) const
      {
        return hostEntityIterator_ == rhs.hostEntityIterator_;
      }

      Entity dereference () const
      {
        return EntityImpl( grid(), *hostEntityIterator_ );
      }

      int level () const { return hostEntityIterator_.level(); }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

    private:
      const Grid *grid_;
      HostEntityIterator hostEntityIterator_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ITERATOR_HH
