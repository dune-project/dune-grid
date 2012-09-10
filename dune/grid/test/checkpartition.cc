// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKPARTITION_CC
#define DUNE_GRID_TEST_CHECKPARTITION_CC

#include <map>

#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/entitycommhelper.hh>

template< Dune::PartitionIteratorType pitype >
struct PartitionFilter;

template<>
struct PartitionFilter< Dune::Interior_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::InteriorEntity);
  }
};

template<>
struct PartitionFilter< Dune::InteriorBorder_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::InteriorEntity) || (partitionType == Dune::BorderEntity);
  }
};

template<>
struct PartitionFilter< Dune::Overlap_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType != Dune::FrontEntity) && (partitionType != Dune::GhostEntity);
  }
};

template<>
struct PartitionFilter< Dune::OverlapFront_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType != Dune::GhostEntity);
  }
};

template<>
struct PartitionFilter< Dune::All_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return true;
  }
};

template<>
struct PartitionFilter< Dune::Ghost_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::GhostEntity);
  }
};



inline bool possibleSubPartitionType ( Dune::PartitionType ept, Dune::PartitionType pt )
{
  switch( ept )
  {
  case Dune::InteriorEntity :
    return (pt == Dune::InteriorEntity) || (pt == Dune::BorderEntity);
  case Dune::OverlapEntity :
    return (pt == Dune::BorderEntity) || (pt == Dune::OverlapEntity) || (pt == Dune::FrontEntity);
  case Dune::GhostEntity :
    return (pt == Dune::BorderEntity) || (pt == Dune::FrontEntity) || (pt == Dune::GhostEntity);
  default :
    std::cerr << "Error: Codimension 0 entity cannot be of partition type " << ept << "." << std::endl;
    return false;
  }
}



template< class GridView, Dune::PartitionIteratorType pitype >
class CheckPartitionType
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking iterators for " << pitype << "..." << std::endl;
    Dune::ForLoop< CheckCodim, 0, GridView::dimension >::apply( gridView );
  }
};


template< class GridView, Dune::PartitionIteratorType pitype >
template< int codim >
struct CheckPartitionType< GridView, pitype >::CheckCodim
{
  typedef typename GridView::template Codim< codim >::template Partition< pitype >::Iterator Iterator;
  typedef typename GridView::template Codim< 0 >::template Partition< Dune::All_Partition >::Iterator AllIterator;

  template< class IdSet >
  static void check ( const Dune::true_type &, const GridView &gridView,
                      const IdSet &idSet )
  {
    typedef std::map< typename IdSet::IdType, Dune::PartitionType > Map;
    typedef typename Map::iterator MapIterator;
    Map map;

    try {
      const Iterator end = gridView.template end< codim, pitype >();
      for( Iterator it = gridView.template begin< codim, pitype >(); it != end; ++it )
      {
        Dune::PartitionType pt = it->partitionType();
        if( !PartitionFilter< pitype >::contains( pt ) )
        {
          std::cerr << "Error: Codim " << codim << " iterator for the " << pitype
                    << " visited entity " << idSet.id( *it )
                    << " with partition type " << pt << "." << std::endl;
        }
        map[ idSet.id( *it ) ] = pt;
      }

      const AllIterator allEnd = gridView.template end< 0, Dune::All_Partition >();
      for( AllIterator it = gridView.template begin< 0, Dune::All_Partition >(); it != allEnd; ++it )
      {
        Dune::PartitionType ept = it->partitionType();
        for( int i = 0; i < it->template count< codim >(); ++i )
        {
          Dune::PartitionType pt = it->template subEntity< codim >( i )->partitionType();
          if( !possibleSubPartitionType( ept, pt ) )
          {
            std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                      << " with partition type " << pt << " is a subentity of entity "
                      << idSet.id( *it ) << " with partition type " << ept << "." << std::endl;
          }

          const MapIterator mapIt = map.find( idSet.subId( *it, i, codim ) );
          if( mapIt == map.end() )
          {
            if( PartitionFilter< pitype >::contains( pt ) )
            {
              std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                        << " with partition type " << pt << " is not visited by codim " << codim
                        << " iterator for the " << pitype << "." << std::endl;
            }
          }
          else
          {
            if( pt != mapIt->second )
            {
              std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                        << " with partition type " << pt << " reported partition type "
                        << mapIt->second << " when accessed with the codim " << codim
                        << " iterator for the " << pitype << "." << std::endl;
            }
          }
        }
      }
    }
    catch( const Dune::Exception &exception )
    {
      std::cerr << "Error: Caught  exception when testing partition iterator for the " << pitype
                << "(" << exception << ")." << std::endl;
    }
  }

  template< class IdSet >
  static void check ( const Dune::false_type &, const GridView &gridView, const IdSet &idSet )
  {}

  static void apply ( const GridView &gridView )
  {
    Dune::integral_constant<
        bool, Dune::Capabilities::hasEntity< typename GridView::Grid, codim >::v
        > capabilityVariable;
    check( capabilityVariable, gridView, gridView.grid().localIdSet() );
  }
};



template< class GridView, Dune::InterfaceType iftype  >
class CheckPartitionDataHandle
  : public Dune::CommDataHandleIF< CheckPartitionDataHandle< GridView, iftype >, int >
{
  static const int dimension = GridView::dimension;

  typedef typename GridView::Grid Grid;

  template< int codim >
  struct Contains
  {
    static void apply ( bool (&contains)[ dimension+1 ] )
    {
      contains[ codim ] = Dune::Capabilities::canCommunicate< Grid, codim >::v;
    }
  };

public:
  explicit CheckPartitionDataHandle ( const Grid &grid )
    : grid_( grid ),
      rank_( grid_.comm().rank() ),
      invalidDimension_( false ),
      invalidCodimension_( false ),
      invalidEntity_( false ),
      invalidSendEntity_( false ),
      invalidReceiveEntity_( false ),
      invalidSize_( false ),
      selfReceive_( false ),
      doubleInterior_( false )
  {
    Dune::ForLoop< Contains, 0, dimension >::apply( contains_ );
  }

  ~CheckPartitionDataHandle ()
  {
    if( invalidDimension_ )
      std::cerr << "[ " << rank_ << " ] Error: Invalid dimension passed during communication." << std::endl;
    if( invalidCodimension_ )
      std::cerr << "[ " << rank_ << " ] Error: Invalid codimension passed during communication." << std::endl;
    if( invalidEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Uncontained entity passed during communication." << std::endl;
    if( invalidSendEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Sent data on entity not contained in communication interface." << std::endl;
    if( invalidReceiveEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Received data on entity not contained in communication interface." << std::endl;
    if( invalidSize_ )
      std::cerr << "[ " << rank_ << " ] Error: Wrong size passed during communication." << std::endl;
    if( selfReceive_ )
      std::cerr << "[ " << rank_ << " ] Warning: Received data from own process during communication." << std::endl;
    if( doubleInterior_ )
      std::cerr << "[ " << rank_ << " ] Error: Received interior data on interior entity." << std::endl;
  }

  bool contains ( int const dim, const int codim ) const
  {
    invalidDimension_ |= (dim != dimension);
    invalidCodimension_ |= ((codim < 0) || (codim > dimension));
    return ((codim >= 0) && (codim <= dimension) ? contains_[ codim ] : false);
  }

  bool fixedsize ( const int dim, const int codim ) const
  {
    invalidDimension_ |= (dim != dimension);
    invalidCodimension_ |= ((codim < 0) || (codim > dimension));
    return true;
  }

  template< class Entity >
  size_t size ( const Entity &entity ) const
  {
    dune_static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    dune_static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );
    return (contains_[ Entity::codimension ] ? 2 : 0);
  }

  template< class Buffer, class Entity >
  void gather ( Buffer &buffer, const Entity &entity ) const
  {
    dune_static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    dune_static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );

    invalidEntity_ |= !contains_[ Entity::codimension ];
    invalidSendEntity_ |= !Dune::EntityCommHelper< iftype >::send( entity.partitionType() );

    buffer.write( rank_ );
    buffer.write( int( entity.partitionType() ) );
  }

  template< class Buffer, class Entity >
  void scatter ( Buffer &buffer, const Entity &entity, size_t n )
  {
    dune_static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    dune_static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );

    invalidEntity_ |= !contains_[ Entity::codimension ];
    invalidSize_ |= (n != size( entity ));
    invalidSendEntity_ |= !Dune::EntityCommHelper< iftype >::receive( entity.partitionType() );

    int rank, partitionType;
    buffer.read( rank );
    buffer.read( partitionType );

    selfReceive_ |= (rank == rank_);
    if( (partitionType == int( Dune::InteriorEntity )) && (entity.partitionType() == Dune::InteriorEntity) )
    {
      std::cout << "[ " << rank_ << " ] Error: Receive interior data from process " << rank
                << " on interior entity " << grid_.globalIdSet().id( entity ) << "." << std::endl;
      doubleInterior_ = true;
    }
  }

  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking communication for " << iftype << "..." << std::endl;
    CheckPartitionDataHandle handle( gridView.grid() );
    gridView.communicate( handle, iftype, Dune::ForwardCommunication );
  }

private:
  const Grid &grid_;
  const int rank_;
  bool contains_[ dimension+1 ];
  mutable bool invalidDimension_;
  mutable bool invalidCodimension_;
  mutable bool invalidEntity_;
  mutable bool invalidSendEntity_;
  bool invalidReceiveEntity_;
  bool invalidSize_;
  bool selfReceive_;
  bool doubleInterior_;
};



template< class GridView >
inline void checkPartitionType ( const GridView &gridView )
{
  CheckPartitionType< GridView, Dune::Interior_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::InteriorBorder_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::Overlap_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::OverlapFront_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::Ghost_Partition >::apply( gridView );

  CheckPartitionDataHandle< GridView, Dune::InteriorBorder_InteriorBorder_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::InteriorBorder_All_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::Overlap_OverlapFront_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::Overlap_All_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::All_All_Interface >::apply( gridView );
}

#endif // DUNE_GRID_TEST_CHECKPARTITION_CC
