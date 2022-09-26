// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKDGF_HH
#define DUNE_GRID_TEST_CHECKDGF_HH

#define CHECK 1

#ifndef DUNE_GRID_EXAMPLE_GRIDS_PATH
#define DUNE_GRID_EXAMPLE_GRIDS_PATH
#endif

#include <string>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/grid/test/checkintersectionit.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune
{

  template< int dim, int dimworld >
  class AlbertaGrid;

}

template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};

template< class GridView >
void display ( const std::string &name,
               const GridView &view,
               std::vector< double > &elDat, int nofElParams,
               std::vector< double > &vtxDat, int nofVtxParams )
{
  Dune::VTKWriter<GridView> vtkWriter(view);
  // SubsamplingVTKWriter<GridView> vtkWriter(view,6);
  if( nofElParams + nofVtxParams > 0 )
  {
    vtkWriter.addCellData( elDat, "el. Parameters", nofElParams );
    vtkWriter.addVertexData( vtxDat, "vtx. Parameters", nofVtxParams );
  }
  auto pos = name.find_last_of("\\/") + 1;
  vtkWriter.write( name.substr(pos, name.size() - pos) );
}

template< class Grid >
void test ( Grid &grid )
{
  gridcheck( grid );

  // check the method geometryInFather()
  std::cout << "  CHECKING: geometry in father" << std::endl;
  checkGeometryInFather( grid );
  // check the intersection iterator and the geometries it returns
  std::cout << "  CHECKING: intersections" << std::endl;
  bool skip = ! EnableLevelIntersectionIteratorCheck< Grid >::v;
  checkIntersectionIterator( grid, skip );
}

template<typename GridType>
void runDGFTest(int argc, char ** argv)
{
  using namespace Dune;

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
#ifdef DGFTEST_USE_GMSH
    namestr << "hybrid-testgrid-" << GridType::dimension << "d.msh";
    filename = std::string( DUNE_GRID_EXAMPLE_GRIDS_PATH ) + "gmsh/" + namestr.str();
#else
    namestr << "test" << GridType::dimension << "d.dgf";
    filename = std::string( DUNE_GRID_EXAMPLE_GRIDS_PATH ) + "dgf/" + namestr.str();
#endif
  }

  std::cout << "tester: start grid reading; file " << filename << std::endl;

  typedef typename GridType::LeafGridView GridView;
  typedef typename GridView::IndexSet IndexSetType;

  // create Grid from DGF parser
  GridType *grid;
  size_t nofElParams( 0 ), nofVtxParams( 0 );
  std::vector< double > eldat( 0 ), vtxdat( 0 );
  {
    GridPtr< GridType > gridPtr( filename, mpiHelper.getCommunicator() );

    gridPtr.loadBalance();

    GridView gridView = gridPtr->leafGridView();
    const IndexSetType &indexSet = gridView.indexSet();
    nofElParams = gridPtr.nofParameters( 0 );
    if( nofElParams > 0 )
    {
      std::cout << "Reading Element Parameters:" << std::endl;
      eldat.resize( indexSet.size(0) * nofElParams );
      const PartitionIteratorType partType = All_Partition;
      typedef typename GridView::template Codim< 0 >::template Partition< partType >::Iterator Iterator;
      const Iterator enditer = gridView.template end< 0, partType >();
      for( Iterator iter = gridView.template begin< 0, partType >(); iter != enditer; ++iter )
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

    nofVtxParams = gridPtr.nofParameters( (int) GridType::dimension );
    if( nofVtxParams > 0 )
    {
      std::cout << "Reading Vertex Parameters:" << std::endl;
      vtxdat.resize( indexSet.size( GridType::dimension ) * nofVtxParams );
      for(const auto & e : elements(gridView /*, Partitions::All */))
      {
        const std::vector< double > &param = gridPtr.parameters( e );
        assert( param.size() == nofVtxParams );
        for( size_t i = 0; i < nofVtxParams; ++i )
        {
          vtxdat[ indexSet.index(e) * nofVtxParams + i ] = param[ i ];
        }
      }
    }


#ifndef ModelP  // parallel UGGrid can only have one grid at a time. defined ModelP => parallel UG
    // does a second construction of the grid work?
    GridPtr< GridType > gridPtr1( filename.c_str(), mpiHelper.getCommunicator() );
#endif

    grid = gridPtr.release();
  }

  GridView gridView = grid->leafGridView();
  std::cout << "Grid size: " << grid->size(0) << std::endl;
  // display
  display( filename , gridView, eldat, nofElParams, vtxdat, nofVtxParams );
  // refine
#if CHECK
  std::cout << "tester: refine grid" << std::endl;
  grid->globalRefine(Dune::DGFGridInfo<GridType>::refineStepsForHalf());
  std::cout << "Grid size: " << grid->size(0) << std::endl;
  test(*grid);
#endif
  delete grid;
}

#endif // #ifndef DUNE_GRID_TEST_CHECKDGF_HH
