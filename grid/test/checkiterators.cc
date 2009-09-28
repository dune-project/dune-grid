// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <map>

#include <dune/grid/genericgeometry/forloop.hh>

template< class GridView >
class CheckIterators
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking iterators for higher codimension..." << std::endl;
    ForLoop< CheckCodim, 1, GridView::dimension >::apply( gridView );
  }
};



template< class GridView >
template< int codim >
struct CheckIterators< GridView >::CheckCodim
{
  typedef typename GridView::Grid::Traits::LocalIdSet LocalIdSet;
  typedef typename LocalIdSet::IdType IdType;

  typedef typename GridView::template Codim< codim >::Iterator CodimIterator;
  typedef typename GridView::template Codim< 0 >::Iterator ElementIterator;

  static void apply ( const GridView &gridView )
  {
    const LocalIdSet &idSet = gridView.grid().localIdSet();

    std::map< IdType, int > count;
    int size = 0;

    const CodimIterator codimEnd = gridView.template end< codim >();
    for( CodimIterator it = gridView.template begin< codim >(); it != codimEnd; ++it )
    {
      ++count[ idSet.id( *it ) ];
      ++size;
    }

    const ElementIterator elementEnd = gridView.template end< 0 >();
    for( ElementIterator it = gridView.template begin< 0 >(); it != elementEnd; ++it )
    {
      const typename ElementIterator::Entity &entity = *it;
      for( int i = 0; i < entity.template count< codim >(); ++i )
      {
        IdType idOld = idSet.template subId<codim>( entity, i);
        idOld = 0;
        IdType id = idSet.subId( entity, i, codim );
        if( count[ id ] != 1 )
        {
          std::cerr << "Error: Codim " << codim << " iterator"
                    << " visited entity " << id
                    << " " << count[ id ] << " times." << std::endl;
          assert( count[ id ] == 1 );
        }
      }
    }
  }
};



template< class GridView >
inline void checkIterators ( const GridView &gridView )
{
  CheckIterators< GridView >::apply( gridView );
}
