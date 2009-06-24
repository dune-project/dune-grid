// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/smartpointer.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#endif



// TriangulatedUnitSquare

template< typename Grid >
class TriangulatedUnitSquareMaker
{
  dune_static_assert( Grid::dimension == 2, "Dimension of grid must be 2" );
  dune_static_assert( Grid::dimensionworld == 2, "Dimension of world must be 2" );

public:
  static Dune::SmartPointer< Grid > create ()
  {
    Dune::GridFactory< Grid > gf;
    Dune::FieldVector< typename Grid::ctype, 2 > pos;

    pos[ 0 ] = 0; pos[ 1 ] = 0; gf.insertVertex( pos );
    pos[ 0 ] = 1; pos[ 1 ] = 0; gf.insertVertex( pos );
    pos[ 0 ] = 0; pos[ 1 ] = 1; gf.insertVertex( pos );
    pos[ 0 ] = 1; pos[ 1 ] = 1; gf.insertVertex( pos );

    Dune::GeometryType type;
    type.makeTriangle();
    std::vector< unsigned int > vid( 3 );

    vid[ 0 ] = 0; vid[ 1 ] = 1; vid[ 2 ] = 2; gf.insertElement( type, vid );
    vid[ 0 ] = 1; vid[ 1 ] = 2; vid[ 2 ] = 3; gf.insertElement( type, vid );

    return gf.createGrid();
  }
};

#ifdef HAVE_ALUGRID
template<>
class TriangulatedUnitSquareMaker< Dune::ALUSimplexGrid< 2, 2 > >
{
  typedef Dune::ALUSimplexGrid< 2, 2 > Grid;

public:
  static Dune::SmartPointer< Grid > create ()
  {
    return new Grid( "2dsimplex.alu" );
  }
};
#endif // HAVE_ALUGRID



// test

template< typename Grid >
void test( Dune::SmartPointer< Grid > grid, int &result, const std::string &name, unsigned int refine = 0 )
{
  typedef typename Grid::LeafGridView GridView;
  typedef typename GridView::IndexSet IndexSet;
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;

  grid->globalRefine( refine );
  GridView gridView = grid->leafView();
  const IndexSet &indexSet = gridView.indexSet();

  const unsigned int numEdges = indexSet.size( 1 );
  const unsigned int numElements = indexSet.size( 0 );
  std::vector< std::vector< double > > data_( numEdges );
  for( unsigned int i = 0; i < numEdges; ++i )
    data_[ i ].resize( numElements, double( 0 ) );

  const Iterator end = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    const typename Iterator::Entity &entity = *it;
    const unsigned int index = indexSet.index( entity );
    for( int k = 0; k < entity.template count< 1 >(); ++k )
    {
      //data_[ indexSet.template subIndex< 1 >( entity, k ) ][ index ] = double( 1 );
      data_[ indexSet.subIndex( entity, k, 1 ) ][ index ] = double( 1 );
    }
  }

  Dune::VTKWriter< GridView > vtkWriter( gridView );
  for( unsigned int i = 0; i < numEdges; ++i )
  {
    std::ostringstream name;
    name << "edge_" << i;
    vtkWriter.addCellData( data_[ i ], name.str() );
  }
  std::ostringstream filename;
  filename << "2d-edges-" << name;
  vtkWriter.write( filename.str(), Dune::VTKOptions::ascii );

  result = 0;
}



// main

int main( int argc, char** argv )
try
{
  //Maybe initialize Mpi
  Dune::MPIHelper::instance( argc, argv );

  // default exitcode 77 (=skipped);
  // returned in case none of the supported Grids were found
  int result = 77;

#ifdef HAVE_ALBERTA
  test( TriangulatedUnitSquareMaker< Dune::AlbertaGrid< 2, 2 > >::create(),
        result, "alberta", 2 );
#endif

#ifdef HAVE_ALUGRID
  test( TriangulatedUnitSquareMaker< Dune::ALUSimplexGrid<2, 2> >::create(),
        result, "alu", 1 );
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
  test( TriangulatedUnitSquareMaker< Dune::UGGrid< 2 > >::create(),
        result, "ug", 1 );
#endif // HAVE_UG

  return result;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
