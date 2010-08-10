// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define DISABLE_DEPRECATED_METHOD_CHECK 1
#define NEW_SUBENTITY_NUMBERING 1
#define CHECK 1

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

/*
   namespace Dune
   {

   template< int dim, int dimworld >
   class AlbertaGrid;

   }

   template< int dim, int dimworld >
   struct EnableLevelIntersectionIteratorCheck< AlbertaGrid< dim, dimworld > >
   {
   static const bool v = false;
   };
 */
using namespace Dune;
template< class GridView >
void display ( const std::string &name,
               const GridView &view,
               std::vector< double > &elDat, int nofElParams,
               std::vector< double > &vtxDat, int nofVtxParams )
{
#if HAVE_GRAPE
  if( nofElParams + nofVtxParams > 0 )
  {
    GrapeDataDisplay< typename GridView::Grid > disp( view );
    disp.addVector( "el. Paramters", elDat, view.indexSet(), 0.0, 0, nofElParams, false );
    disp.addVector( "vtx. Paramters", vtxDat, view.indexSet(), 0.0, 1, nofVtxParams, true );
    disp.display();
  }
  else
  {
    GrapeGridDisplay< typename GridView::Grid > disp( view.grid() );
    disp.display();
  }
#endif // #if HAVE_GRAPE
  VTKWriter<GridView> vtkWriter(view);
  // SubsamplingVTKWriter<GridView> vtkWriter(view,6);
  if( nofElParams + nofVtxParams > 0 )
  {
    vtkWriter.addCellData( elDat, "el. Parameters", nofElParams );
    vtkWriter.addVertexData( vtxDat, "vtx. Parameters", nofVtxParams );
  }
  vtkWriter.write( name );
}

template< class GridView >
void test ( const GridView &view )
{
  typename GridView::Grid &grid = const_cast< typename GridView::Grid & >( view.grid());
  gridcheck( grid );

  // check the method geometryInFather()
  std::cout << "  CHECKING: geometry in father" << std::endl;
  checkGeometryInFather( grid );
  // check the intersection iterator and the geometries it returns
  std::cout << "  CHECKING: intersections" << std::endl;
  checkIntersectionIterator( grid );
}

int main(int argc, char ** argv, char ** envp)
try {
  typedef GridSelector::GridType GridType;
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);

  std::string filename;

  if (argc>=2)
  {
    filename = argv[1];
  }
  else
  {
    std::stringstream namestr;
    namestr << "test" << GridType::dimension << "d.dgf";
    filename = std::string( DUNE_GRID_EXAMPLE_GRIDS_PATH ) + "dgf/" + namestr.str();
  }

  std::cout << "tester: start grid reading; file " << filename << std::endl;

  typedef GridType::LeafGridView GridView;
  typedef GridView::IndexSet IndexSetType;

  // create Grid from DGF parser
  GridType *grid;
  size_t nofElParams( 0 ), nofVtxParams( 0 );
  std::vector< double > eldat( 0 ), vtxdat( 0 );
  {
    GridPtr< GridType > gridPtr( filename.c_str(), mpiHelper.getCommunicator() );

    gridPtr.loadBalance();

    GridView gridView = gridPtr->leafView();
    const IndexSetType &indexSet = gridView.indexSet();
    nofElParams = gridPtr.nofParameters( 0 );
    if( nofElParams > 0 )
    {
      std::cout << "Reading Element Parameters:" << std::endl;
      eldat.resize( indexSet.size(0) * nofElParams );
      const PartitionIteratorType partType = All_Partition;
      typedef GridView::Codim< 0 >::Partition< partType >::Iterator Iterator;
      const Iterator enditer = gridView.end< 0, partType >();
      for( Iterator iter = gridView.begin< 0, partType >(); iter != enditer; ++iter )
      {
        const std::vector< double > &param = gridPtr.parameters( *iter );
        assert( param.size() == nofElParams );
        for( size_t i = 0; i < nofElParams; ++i )
        {
          //std::cout << param[i] << " ";
          eldat[ indexSet.index(*iter) * nofElParams + i ] = param[ i ];
        }
        //std::cout << std::endl;
      }
    }

    nofVtxParams = gridPtr.nofParameters( GridType::dimension );
    if( nofVtxParams > 0 )
    {
      std::cout << "Reading Vertex Parameters:" << std::endl;
      vtxdat.resize( indexSet.size( GridType::dimension ) * nofVtxParams );
      const PartitionIteratorType partType = All_Partition;
      typedef GridView::Codim< GridType::dimension >::Partition< partType >::Iterator Iterator;
      const Iterator enditer = gridView.end< GridType::dimension, partType >();
      for( Iterator iter = gridView.begin< GridType::dimension, partType >(); iter != enditer; ++iter )
      {
        const std::vector< double > &param = gridPtr.parameters( *iter );
        assert( param.size() == nofVtxParams );
        // std::cout << (*iter).geometry()[0] << " -\t ";
        for( size_t i = 0; i < nofVtxParams; ++i )
        {
          // std::cout << param[i] << " ";
          vtxdat[ indexSet.index(*iter) * nofVtxParams + i ] = param[ i ];
        }
        // std::cout << std::endl;
      }
    }

    // does a second construction of the grid work?
    GridPtr< GridType > gridPtr1( filename.c_str(), mpiHelper.getCommunicator() );

    grid = gridPtr.release();
  }

  GridView gridView = grid->leafView();
  std::cout << "Grid size: " << grid->size(0) << std::endl;
  // display
  display( filename , gridView, eldat, nofElParams, vtxdat, nofVtxParams );
  // refine
#if CHECK
  std::cout << "tester: refine grid" << std::endl;
  grid->globalRefine(Dune::DGFGridInfo<GridType>::refineStepsForHalf());
  std::cout << "Grid size: " << grid->size(0) << std::endl;
  test(gridView);
#endif
  delete grid;
  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
