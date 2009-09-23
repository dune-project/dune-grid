// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <map>

#include <dune/grid/genericgeometry/misc.hh>

using GenericGeometry::ForLoop;


template< PartitionIteratorType pitype >
struct PartitionFilter;

template<>
struct PartitionFilter< Interior_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return (partitionType == InteriorEntity);
  }
};

template<>
struct PartitionFilter< InteriorBorder_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return (partitionType == InteriorEntity) || (partitionType == BorderEntity);
  }
};

template<>
struct PartitionFilter< Overlap_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return (partitionType == OverlapEntity);
  }
};

template<>
struct PartitionFilter< OverlapFront_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return (partitionType == OverlapEntity) || (partitionType == FrontEntity);
  }
};

template<>
struct PartitionFilter< All_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return true;
  }
};

template<>
struct PartitionFilter< Ghost_Partition >
{
  static bool contains ( const PartitionType partitionType )
  {
    return (partitionType == GhostEntity);
  }
};



inline bool possibleSubPartitionType ( PartitionType ept, PartitionType pt )
{
  switch( ept )
  {
  case InteriorEntity :
    return (pt == InteriorEntity) || (pt == BorderEntity);
  case OverlapEntity :
    return (pt == BorderEntity) || (pt == OverlapEntity) || (pt == FrontEntity);
  case GhostEntity :
    return (pt == BorderEntity) || (pt == FrontEntity) || (pt == GhostEntity);
  default :
    std::cerr << "Error: Codimension 0 entity cannot be of partition type " << ept << "." << std::endl;
    return false;
  }
}



template< class GridView, PartitionIteratorType pitype >
class CheckPartitionType
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking partition iterators..." << std::endl;
    ForLoop< CheckCodim, 0, GridView::dimension >::apply( gridView );
  }
};


template< class GridView, PartitionIteratorType pitype >
template< int codim >
struct CheckPartitionType< GridView, pitype >::CheckCodim
{
  typedef typename GridView::Grid::Traits::LocalIdSet LocalIdSet;
  typedef typename LocalIdSet::IdType IdType;

  typedef typename GridView::template Codim< codim >::template Partition< pitype >::Iterator Iterator;
  typedef typename GridView::template Codim< 0 >::template Partition< All_Partition >::Iterator AllIterator;

  static void apply ( const GridView &gridView )
  {
    const LocalIdSet &idSet = gridView.grid().localIdSet();

    typedef std::map< IdType, PartitionType > Map;
    typedef typename Map::iterator MapIterator;
    Map map;

    const Iterator end = gridView.template end< codim, pitype >();
    for( Iterator it = gridView.template begin< codim, pitype >(); it != end; ++it )
    {
      PartitionType pt = it->partitionType();
      if( !PartitionFilter< pitype >::contains( pt ) )
      {
        std::cerr << "Error: Codim " << codim << " iterator for the " << pitype
                  << " visited entity " << idSet.id( *it )
                  << " with partition type " << pt << "." << std::endl;
      }
      map[ idSet.id( *it ) ] = pt;
    }

    const AllIterator allEnd = gridView.template end< 0, All_Partition >();
    for( AllIterator it = gridView.template begin< 0, All_Partition >(); it != allEnd; ++it )
    {
      PartitionType ept = it->partitionType();
      for( int i = 0; i < it->template count< codim >(); ++i )
      {
        PartitionType pt = it->template subEntity< codim >( i )->partitionType();
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
};



template< class GridView >
inline void checkPartitionType ( const GridView &gridView )
{
  CheckPartitionType< GridView, Interior_Partition >::apply( gridView );
  CheckPartitionType< GridView, InteriorBorder_Partition >::apply( gridView );
  CheckPartitionType< GridView, Overlap_Partition >::apply( gridView );
  CheckPartitionType< GridView, OverlapFront_Partition >::apply( gridView );
  CheckPartitionType< GridView, Ghost_Partition >::apply( gridView );
}
