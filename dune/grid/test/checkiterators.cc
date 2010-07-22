// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <map>

#include <dune/common/forloop.hh>


// CheckCodimIterators
// -------------------

template< class GridView, int codim,
    bool hasEntity = Dune::Capabilities::hasEntity< typename GridView::Grid, codim >::v >
struct CheckCodimIterators;

template< class GridView, int codim >
struct CheckCodimIterators< GridView, codim, false >
{
  static void apply ( const GridView &gridView )
  {
    std::cerr << "Warning: Not checking iterators for codimension " << codim
              << ", because the corresponding entities are not implemented." << std::endl;
  }
};

template< class GridView, int codim >
struct CheckCodimIterators< GridView, codim, true >
{
  static void apply ( const GridView &gridView );
};


// CheckIterators
// --------------

template< class GridView >
class CheckIterators
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking iterators for higher codimension..." << std::endl;
    Dune::ForLoop< CheckCodim, 1, GridView::dimension >::apply( gridView );
  }
};


// CheckIterators::CheckCodim
// --------------------------

template< class GridView >
template< int codim >
struct CheckIterators< GridView >::CheckCodim
{
  static void apply ( const GridView &gridView )
  {
    return CheckCodimIterators< GridView, codim >::apply( gridView );
  }
};


// Implementation of CheckCodimIterators
// -------------------------------------

template< class GridView, int codim >
inline void CheckCodimIterators< GridView, codim, true >
::apply ( const GridView &gridView )
{
  typedef typename GridView::Grid::Traits::LocalIdSet LocalIdSet;
  typedef typename LocalIdSet::IdType IdType;

  typedef typename GridView::template Codim< codim >::Iterator CodimIterator;
  typedef typename GridView::template Codim< 0 >::Iterator ElementIterator;

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



template< class GridView >
inline void checkIterators ( const GridView &gridView )
{
  CheckIterators< GridView >::apply( gridView );
}
