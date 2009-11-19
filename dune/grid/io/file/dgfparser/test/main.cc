// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define DISABLE_DEPRECATED_METHOD_CHECK 1
#define NEW_SUBENTITY_NUMBERING 1

#include <dune/grid/test/gridcheck.cc>
#include "../dgfgridtype.hh"

namespace Dune
{

  template< int dim, int dimworld >
  class AlbertaGrid;

}

using namespace Dune;
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>

template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};

template< class GridView >
void test ( const GridView &view,
            std::vector< double > &elDat, int nofElParams,
            std::vector< double > &vtxDat, int nofVtxParams )
{
  gridcheck( const_cast< typename GridView::Grid & >( view.grid() ) );
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
}
#else
template< class GridView >
void test( const GridView &view,std::vector<double>& dat,int nofParams,int cdim)
{
  gridcheck( const_cast< typename GridView::Grid & >( view.grid() ) );
}
#endif
template< class GridView >
void test ( const GridView &view )
{
  gridcheck( const_cast< typename GridView::Grid & >( view.grid() ) );
}

int main(int argc, char ** argv, char ** envp)
try {
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);
  int myrank = mpiHelper.rank();

  if (argc<2) {
    std::cerr << "supply grid file as parameter!" << std::endl;
    return 1;
  }

  std::cout << "tester: start grid reading" << std::endl;

  // create Grid from DGF parser
  GridPtr<GridType> grid;
  {
    grid = GridPtr<GridType>(argv[1], mpiHelper.getCommunicator() );
  }

  typedef GridType::LeafGridView GridView;
  typedef GridView::IndexSet IndexSetType;

  GridView gridView = grid->leafView();
  const IndexSetType &index = gridView.indexSet();
  std::vector<double> eldat(0),vtxdat(0);
  if (grid.nofParameters(0)>0)
  {
    std::cout << "Reading Element Parameters:" << std::endl;
    eldat.resize(index.size(0)*grid.nofParameters(0));
    typedef GridView::Codim< 0 >::Iterator IteratorType;
    const IteratorType endit = gridView.end< 0 >();
    for (IteratorType iter=gridView.begin< 0 >(); iter!=endit; ++iter)
    {
      std::vector<double>& param = grid.parameters(*iter);
      assert( (int) param.size() == grid.nofParameters(0) );
      for (size_t i=0; i<param.size(); i++)
      {
        //std::cout << param[i] << " ";
        eldat[index.index(*iter)*param.size()+i] = param[i];
      }
      //std::cout << std::endl;
    }
  }
  if (grid.nofParameters(GridType::dimension)>0)
  {
    std::cout << "Reading Vertex Parameters:" << std::endl;
    vtxdat.resize(index.size(GridType::dimension)*grid.nofParameters(GridType::dimension));
    typedef GridView::Codim< GridType::dimension >::Iterator IteratorType;
    IteratorType endit = gridView.end< GridType::dimension >();
    for (IteratorType iter=gridView.begin< GridType::dimension >();
         iter!=endit; ++iter) {
      std::vector<double>& param = grid.parameters(*iter);
      assert( (int) param.size() == grid.nofParameters(GridType::dimension));
      // std::cout << (*iter).geometry()[0] << " -\t ";
      for (size_t i=0; i<param.size(); i++) {
        // std::cout << param[i] << " ";
        vtxdat[index.index(*iter)*param.size()+i] = param[i];
      }
      // std::cout << std::endl;
    }
  }

  // display
  if (myrank <= 0)
    test(gridView,eldat,grid.nofParameters(0),vtxdat,grid.nofParameters(GridType::dimension));
  // refine
  std::cout << "tester: refine grid" << std::endl;
  grid->globalRefine(Dune::DGFGridInfo<GridType>::refineStepsForHalf());
  test(gridView);
  return 0;
}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
}
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
