// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

#include <dune/grid/common/referenceelements.hh>

#include <dune/grid/geogrid/entitypointer.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< class Traits, bool fake = Traits :: fake >
  class GeometryGridIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator;

  template< class Grid >
  class GeometryGridHierarchicIterator;



  // GeometryGridPartitionIteratorFilter
  // -----------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct GeometryGridPartitionIteratorFilter;

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, Interior_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = Interior_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      const int size = refElement.size( subEntity, codim, dimension );
      for( int i = 0; i < size; ++i )
      {
        const int j = refElement.subEntity( subEntity, codim, i, dimension );
        PartitionType type = element.template entity< dimension >( j )->partitionType();
        if( type == InteriorEntity )
          return true;
      }
      return false;
    }
  };

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, InteriorBorder_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = Interior_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      return true;
    }
  };

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, Overlap_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = Overlap_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      const int size = refElement.size( subEntity, codim, dimension );
      bool border = false;
      bool front = false;
      for( int i = 0; i < size; ++i )
      {
        const int j = refElement.subEntity( subEntity, codim, i, dimension );
        PartitionType type = element.template entity< dimension >( j )->partitionType();
        if( type == OverlapEntity )
          return true;
        border |= (type == BorderEntity);
        front |= (type == FrontEntity);
      }
      return (border && front);
    }
  };

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, OverlapFront_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = Overlap_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      const int size = refElement.size( subEntity, codim, dimension );
      for( int i = 0; i < size; ++i )
      {
        const int j = refElement.subEntity( subEntity, codim, i, dimension );
        PartitionType type = element.template entity< dimension >( j )->partitionType();
        if( type != BorderEntity )
          return true;
      }
      return false;
    }
  };

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, All_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = All_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      return true;
    }
  };

  template< int codim, class Grid >
  struct GeometryGridPartitionIteratorFilter< codim, Ghost_Partition, Grid >
  {
    static const int dimension = Grid :: dimension;
    static const int codimension = codim;

    static const PartitionIteratorType Element_Partition = Ghost_Partition;

    typedef typename Grid :: template Codim< 0 > :: Entity Element;
    typedef Dune :: ReferenceElement< typename Grid :: ctype, dimension >
    ReferenceElement;

    static bool apply ( const ReferenceElement &refElement,
                        const Element &element, int subEntity )
    {
      const int size = refElement.size( subEntity, codim, dimension );
      for( int i = 0; i < size; ++i )
      {
        const int j = refElement.subEntity( subEntity, codim, i, dimension );
        PartitionType type = element.template entity< dimension >( j )->partitionType();
        if( type == GhostEntity )
          return true;
      }
      return false;
    }
  };



  // GeometryGridIterator (real)
  // ---------------------------

  template< class Traits >
  class GeometryGridIterator< Traits, false >
    : public GeometryGridEntityPointer< Traits, false >
  {
    typedef GeometryGridEntityPointer< Traits, false > ExactBase;

    typedef typename Traits :: Grid Grid;

  public:
    typedef typename Traits :: IteratorType IteratorType;

  protected:
    using ExactBase :: hostEntityIterator_;
    using ExactBase :: update;

  public:
    GeometryGridIterator ( const Grid &grid, int level, IteratorType type )
      : ExactBase( grid, Traits :: getHostEntityIterator( grid, level, type ) )
    {}

    void increment ()
    {
      ++hostEntityIterator_;
      update();
    }
  };



  // GeometryGridIterator (fake)
  // ---------------------------

  template< class Traits >
  class GeometryGridIterator< Traits, true >
    : public GeometryGridEntityPointer< Traits, true >
  {
    typedef GeometryGridEntityPointer< Traits, true > ExactBase;

    typedef typename Traits :: Grid Grid;

  public:
    static const int dimension = Traits :: dimension;
    static const int codimension = Traits :: codimension;

    typedef typename Traits :: IteratorType IteratorType;

  private:
    typedef typename Traits :: Filter Filter;

    typedef typename Traits :: HostElement HostElement;
    typedef typename Traits :: HostElementIterator HostElementIterator;
    typedef typename Traits :: HostIndexSet HostIndexSet;

    HostElementIterator hostEnd_;
    const HostIndexSet *hostIndexSet_;
    std :: vector< bool > visited_;

  protected:
    using ExactBase :: hostElementIterator_;
    using ExactBase :: subEntity_;
    using ExactBase :: update;

  public:
    GeometryGridIterator ( const Grid &grid, int level, IteratorType type )
      : ExactBase( grid, Traits :: getHostElementIterator( grid, level, type ), -1 ),
        hostEnd_( Traits :: getHostElementIterator( grid, level, Traits :: end ) ),
        hostIndexSet_( &Traits :: getHostIndexSet( grid, level ) )
    {
      if( hostElementIterator_ != hostEnd_ )
      {
        visited_.resize( hostIndexSet_->size( codimension ), false );
        increment();
      }
    }

    void increment ()
    {
      typedef typename Traits :: ctype ctype;

      while( hostElementIterator_ != hostEnd_ )
      {
        const HostElement &hostElement = *hostElementIterator_;

        const ReferenceElement< ctype, dimension > &refElement
          = ReferenceElements< ctype, dimension > :: general( hostElement.type() );

        ++subEntity_;
        const int count = refElement.size( codimension );
        for( ; subEntity_ < count; ++subEntity_ )
        {
          if( !Filter :: apply( refElement, hostElement, subEntity_ ) )
            continue;

          const size_t index
            = hostIndexSet_->template subIndex< codimension >( hostElement, subEntity_ );
          if( !visited_[ index ] )
          {
            visited_[ index ] = true;
            update();
            return;
          }
        }
        ++hostElementIterator_;
        subEntity_ = -1;
      }
      update();
    }
  };



  // GeometryGridLeafIteratorTraits
  // ------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct GeometryGridLeafIteratorTraits
    : public GeometryGridEntityPointerTraits< codim, Grid >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: HostGrid HostGrid;

    typedef GeometryGridPartitionIteratorFilter< codim, pitype, HostGrid > Filter;

    static const PartitionIteratorType Entity_Partition = pitype;
    static const PartitionIteratorType Element_Partition = Filter :: Element_Partition;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< Entity_Partition > :: LeafIterator
    HostEntityIterator;
    typedef typename HostGrid :: template Codim< 0 >
    :: template Partition< Element_Partition > :: LeafIterator
    HostElementIterator;

    typedef typename HostGrid :: Traits :: LeafIndexSet HostIndexSet;

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



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct GeometryGridLeafIterator
    : public GeometryGridIterator
      < GeometryGridLeafIteratorTraits< codim, pitype, Grid > >
  {
    typedef GeometryGridLeafIteratorTraits< codim, pitype, Grid > Traits;

    typedef typename Traits :: IteratorType IteratorType;

    GeometryGridLeafIterator ( const Grid &grid, IteratorType type )
      : GeometryGridIterator< Traits >( grid, -1, type )
    {}
  };



  // GeometryGridLevelIteratorTraits
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct GeometryGridLevelIteratorTraits
    : public GeometryGridEntityPointerTraits< codim, Grid >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: HostGrid HostGrid;

    typedef GeometryGridPartitionIteratorFilter< codim, pitype, HostGrid > Filter;

    static const PartitionIteratorType Entity_Partition = pitype;
    static const PartitionIteratorType Element_Partition = Filter :: Element_Partition;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< Entity_Partition > :: LevelIterator
    HostEntityIterator;
    typedef typename HostGrid :: template Codim< 0 >
    :: template Partition< Element_Partition > :: LevelIterator
    HostElementIterator;

    typedef typename HostGrid :: Traits :: LevelIndexSet HostIndexSet;

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



  // GeometryGridLevelIterator
  // -------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct GeometryGridLevelIterator
    : public GeometryGridIterator
      < GeometryGridLevelIteratorTraits< codim, pitype, Grid > >
  {
    typedef GeometryGridLevelIteratorTraits< codim, pitype, Grid > Traits;

    typedef typename Traits :: IteratorType IteratorType;

    GeometryGridLevelIterator ( const Grid &grid, int level, IteratorType type )
      : GeometryGridIterator< Traits >( grid, level, type )
    {}
  };



  // GeometryGridHierarchicIteratorTraits
  // ------------------------------------

  template< class Grid >
  struct GeometryGridHierarchicIteratorTraits
    : public GeometryGridEntityPointerTraits< 0, Grid >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: HostGrid :: Traits :: HierarchicIterator
    HostEntityIterator;
    typedef typename Traits :: HostGrid :: Traits :: HierarchicIterator
    HostElementIterator;
  };



  // GeometryGridHierarchicIterator
  // ------------------------------

  template< class Grid >
  class GeometryGridHierarchicIterator
    : public GeometryGridEntityPointer< GeometryGridHierarchicIteratorTraits< Grid > >
  {
    typedef GeometryGridHierarchicIteratorTraits< Grid > Traits;

    typedef GeometryGridEntityPointer< Traits > ExactBase;

  protected:
    typedef typename Traits :: HostEntityIterator HostEntityIterator;

    using ExactBase :: hostEntityIterator_;
    using ExactBase :: update;

  public:
    GeometryGridHierarchicIterator ( const Grid &grid,
                                     const HostEntityIterator &hostEntityIterator )
      : ExactBase( grid, hostEntityIterator )
    {}

    void increment ()
    {
      ++hostEntityIterator_;
      update();
    }
  };

}

#endif
