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

  template< int codim, PartitionIteratorType pitype, class Grid,
      bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v) >
  class GeometryGridLeafIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIteratorAdapter;

  template< int codim, PartitionIteratorType pitype, class Grid,
      bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v ) >
  class GeometryGridLevelIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIteratorAdapter;

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



  // GeometryGridLeafIterator (real)
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIterator< codim, pitype, Grid, false >
    : public GeometryGridEntityPointer< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, false > Base;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< pitype > :: LeafIterator
    HostIterator;

    enum IteratorType { begin, end };

  private:
    HostIterator hostIterator_;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLeafIterator ( const Grid &grid, IteratorType type )
      : Base( grid, getHostIterator( grid, type ) ),
        hostIterator_( getHostIterator( grid, type ) )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }

  private:
    static HostIterator getHostIterator ( const Grid &grid, IteratorType type )
    {
      if( type == begin )
        return grid.hostGrid().template leafbegin< codim, pitype >();
      else
        return grid.hostGrid().template leafend< codim, pitype >();
    }
  };



  // GeometryGridLeafIterator (fake)
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIterator< codim, pitype, Grid, true >
    : public GeometryGridEntityPointer< codim, Grid, true >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, true > Base;

    typedef GeometryGridPartitionIteratorFilter< codim, pitype, HostGrid > Filter;

    typedef typename HostGrid :: template Codim< 0 >
    :: template Partition< Filter :: Element_Partition > :: LeafIterator
    HostIterator;
    typedef typename HostGrid :: Traits :: LeafIndexSet HostIndexSet;

    enum IteratorType { begin, end };

  private:
    typedef typename Base :: HostElement HostElement;

    HostIterator hostEndIterator_;
    HostIterator hostIterator_;

    const HostIndexSet *hostIndexSet_;
    std :: vector< bool > visited_;

  protected:
    using Base :: hostElementPointer;
    using Base :: subEntity;
    using Base :: setToTarget;

  public:
    GeometryGridLeafIterator ( const Grid &grid, IteratorType type )
      : Base( grid, getHostIterator( grid, type ), -1 ),
        hostEndIterator_( getHostIterator( grid, end ) ),
        hostIterator_( getHostIterator( grid, type ) ),
        hostIndexSet_( &grid.hostGrid().leafIndexSet() )
    {
      if( hostIterator_ != hostEndIterator_ )
      {
        visited_.resize( hostIndexSet_->size( codim ), false );
        increment();
      }
    }

    void increment ()
    {
      typedef typename Traits :: ctype ctype;
      static const int dimension = Traits :: dimension;

      while( hostIterator_ != hostEndIterator_ )
      {
        const HostElement &hostElement = *hostElementPointer();

        const ReferenceElement< ctype, dimension > &refElement
          = ReferenceElements< ctype, dimension > :: general( hostElement.type() );

        const int count = refElement.size( codim );
        for( int number = subEntity() + 1; number < count; ++number )
        {
          if( !Filter :: apply( refElement, hostElement, number ) )
            continue;

          const size_t index
            = hostIndexSet_->template subIndex< codim >( hostElement, number );
          if( !visited_[ index ] )
          {
            visited_[ index ] = true;
            setToTarget( hostIterator_, number );
            return;
          }
        }
        ++hostIterator_;
        setToTarget( hostIterator_, -1 );
      }
    }

  private:
    static HostIterator getHostIterator ( const Grid &grid, IteratorType type )
    {
      if( type == begin )
        return grid.hostGrid().template leafbegin< 0, pitype >();
      else
        return grid.hostGrid().template leafend< 0, pitype >();
    }
  };



  // GeometryGridLeafIteratorAdapter
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIteratorAdapter
    : public GeometryGridLeafIterator< codim, pitype, Grid >
  {
    typedef GeometryGridLeafIterator< codim, pitype, Grid > BaseType;

  public:
    typedef typename BaseType :: IteratorType IteratorType;

    GeometryGridLeafIteratorAdapter ( const Grid &grid, IteratorType type )
      : BaseType( grid, type )
    {}
  };



  // GeometryGridLeafIterator (real)
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator< codim, pitype, Grid, false >
    : public GeometryGridEntityPointer< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, false > Base;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< pitype > :: LevelIterator
    HostIterator;

    enum IteratorType { begin, end };

  private:
    HostIterator hostIterator_;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLevelIterator ( const Grid &grid, int level, IteratorType type )
      : Base( grid, getHostIterator( grid, level, type ) ),
        hostIterator_( getHostIterator( grid, level, type ) )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }

  private:
    static HostIterator
    getHostIterator ( const Grid &grid, int level, IteratorType type )
    {
      if( type == begin )
        return grid.hostGrid().template lbegin< codim, pitype >( level );
      else
        return grid.hostGrid().template lend< codim, pitype >( level );
    }
  };



  // GeometryGridLevelIterator (fake)
  // --------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator< codim, pitype, Grid, true >
    : public GeometryGridEntityPointer< codim, Grid, true >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, true > Base;

    typedef GeometryGridPartitionIteratorFilter< codim, pitype, HostGrid > Filter;

    typedef typename HostGrid :: template Codim< 0 >
    :: template Partition< Filter :: Element_Partition > :: LevelIterator
    HostIterator;
    typedef typename HostGrid :: Traits :: LevelIndexSet HostIndexSet;

    enum IteratorType { begin, end };

  private:
    typedef typename Base :: HostElement HostElement;

    HostIterator hostEndIterator_;
    HostIterator hostIterator_;

    const HostIndexSet *hostIndexSet_;
    std :: vector< bool > visited_;

  protected:
    using Base :: hostElementPointer;
    using Base :: subEntity;
    using Base :: setToTarget;

  public:
    GeometryGridLevelIterator ( const Grid &grid, int level, IteratorType type )
      : Base( grid, getHostIterator( grid, level, type ), -1 ),
        hostEndIterator_( getHostIterator( grid, level, end ) ),
        hostIterator_( getHostIterator( grid, level, type ) ),
        hostIndexSet_( &grid.hostGrid().levelIndexSet( level ) )
    {
      if( hostIterator_ != hostEndIterator_ )
      {
        visited_.resize( hostIndexSet_->size( codim ), false );
        increment();
      }
    }

    void increment ()
    {
      typedef typename Traits :: ctype ctype;
      static const int dimension = Traits :: dimension;

      while( hostIterator_ != hostEndIterator_ )
      {
        const HostElement &hostElement = *hostElementPointer();

        const ReferenceElement< ctype, dimension > &refElement
          = ReferenceElements< ctype, dimension > :: general( hostElement.type() );

        const int count = refElement.size( codim );
        for( int number = subEntity() + 1; number < count; ++number )
        {
          if( !Filter :: apply( refElement, hostElement, number ) )
            continue;

          const size_t index
            = hostIndexSet_->template subIndex< codim >( hostElement, number );
          if( !visited_[ index ] )
          {
            visited_[ index ] = true;
            setToTarget( hostIterator_, number );
            return;
          }
        }
        ++hostIterator_;
        setToTarget( hostIterator_, -1 );
      }
    }

  private:
    static HostIterator
    getHostIterator ( const Grid &grid, int level, IteratorType type )
    {
      if( type == begin )
        return grid.hostGrid().template lbegin< 0, pitype >( level );
      else
        return grid.hostGrid().template lend< 0, pitype >( level );
    }
  };



  // GeometryGridLevelIteratorAdapter
  // --------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIteratorAdapter
    : public GeometryGridLevelIterator< codim, pitype, Grid >
  {
    typedef GeometryGridLevelIterator< codim, pitype, Grid > BaseType;

  public:
    typedef typename BaseType :: IteratorType IteratorType;

    GeometryGridLevelIteratorAdapter ( const Grid &grid, int level, IteratorType type )
      : BaseType( grid, level, type )
    {}
  };



  // GeometryGridHierarchicIterator
  // ------------------------------

  template< class Grid >
  class GeometryGridHierarchicIterator
    : public GeometryGridEntityPointer< 0, Grid >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< 0, Grid > Base;

    typedef typename HostGrid :: Traits :: HierarchicIterator HostIterator;

  private:
    HostIterator hostIterator_;

  public:
    GeometryGridHierarchicIterator ( const Grid &grid,
                                     const HostIterator &hostIterator )
      : Base( grid, hostIterator ),
        hostIterator_( hostIterator )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };

}

#endif
