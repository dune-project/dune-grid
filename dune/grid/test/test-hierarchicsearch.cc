// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <cstdlib>

#include <array>
#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/grid/yaspgrid.hh>

// UnitCube
// --------

template< class Grid >
class UnitCube;

template< int dimension >
class UnitCube< Dune::YaspGrid< dimension > >
{
public:
  typedef Dune::YaspGrid< dimension > Grid;

  static std::unique_ptr< Grid > create ()
  {
    Dune::FieldVector< double, dimension > domain( 1. );
    std::array< int, dimension > cells;
    cells.fill( 1 );
    return std::make_unique< Grid >( domain, cells );
  }
};



// check
// -----

template< class GridView >
void check ( GridView gridView )
{
  typedef typename GridView::Grid Grid;
  typedef typename GridView::IndexSet IndexSet;

  typedef Dune::HierarchicSearch< Grid, IndexSet > HierarchicSearch;
  HierarchicSearch hsearch( gridView.grid(), gridView.indexSet() );

  typedef typename GridView::template Codim< 0 >::Iterator Iterator;
  typedef typename Iterator::Entity Entity;

  const Iterator end = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    const Entity &entity = *it;
    if( entity != hsearch.findEntity( entity.geometry().center() ) )
      DUNE_THROW( Dune::GridError, "Could not retrieve element in hierarchic search" );
  }
}



int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  // create grid
  const int dimension = 2;
  typedef Dune::YaspGrid< dimension > Grid;
  auto grid = UnitCube< Grid >::create();

  // grid refinement
  int maxLevel = 4;
  if( argc > 1 )
    maxLevel = std::atoi( argv[ 1 ] );
  grid->globalRefine( maxLevel );

  // check hierarchic search
  for( int level = 0; level < grid->maxLevel(); ++level )
    check( grid->levelGridView( level ) );
  check( grid->leafGridView() );

  return 0;
}
catch( Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
catch( std::exception &exception )
{
  std::cerr << exception.what() << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Unknown exception thrown" << std::endl;
  return 1;
}
