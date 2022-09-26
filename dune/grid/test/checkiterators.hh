// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKITERATORS_HH
#define DUNE_GRID_TEST_CHECKITERATORS_HH

#include <iostream>
#include <map>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/iteratortest.hh>

#include <dune/grid/common/capabilities.hh>

template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& ){}
};

// CheckCodimIterators
// -------------------

template< class GridView, int codim,
    bool hasEntity = Dune::Capabilities::hasEntityIterator< typename GridView::Grid, codim >::v >
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
    Dune::Hybrid::forEach( std::make_index_sequence< GridView::dimension >{}, [ & ]( auto i ){ CheckCodim< i+1 >::apply( gridView ); } );
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
  typedef typename GridView::Grid::Traits::LocalIdSet::IdType IdType;

  const auto &idSet = gridView.grid().localIdSet();

  std::map< IdType, int > count;
  int size = 0;

  const auto codimEnd = gridView.template end< codim >();
  for( auto it = gridView.template begin< codim >(); it != codimEnd; ++it )
  {
    ++count[ idSet.id( *it ) ];
    ++size;
  }

  const auto elementEnd = gridView.template end< 0 >();
  for( auto it = gridView.template begin< 0 >(); it != elementEnd; ++it )
  {
    const auto &entity = *it;
    for( std::size_t i = 0; i < entity.subEntities(codim); ++i )
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

  // check forward iterator semantics
  typedef typename GridView::template Codim<codim>::Entity Entity;
  NoopFunctor<Entity> op;
  if(0 != testForwardIterator(gridView.template begin<codim>(),
                               gridView.template end<codim>(), op))
    DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

}



template< class GridView >
inline void checkIterators ( const GridView &gridView )
{
  CheckIterators< GridView >::apply( gridView );
}

#endif // DUNE_GRID_TEST_CHECKITERATORS_HH
