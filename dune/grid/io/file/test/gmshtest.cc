// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

// dune grid includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

// alberta related stuff
#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

#include <dune/grid/onedgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/gmshwriter.hh>

using namespace Dune;

#if HAVE_ALBERTA
template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};
#endif

template <typename GridType>
void testReadingAndWritingGrid( const std::string& path, const std::string& gridName, const std::string& gridManagerName, int refinements)
{
  // Read the grid
  std::cout<<"Using "<<gridManagerName<<std::endl;
  GridFactory<GridType> gridFactory;
  std::vector<int> boundaryIDs;
  std::vector<int> elementsIDs;
  const std::string inputName(path+gridName+".msh");
  std::cout<<"Reading mesh file "<<inputName<<std::endl;
  GmshReader<GridType>::read(gridFactory,inputName,boundaryIDs,elementsIDs);
  auto grid=std::unique_ptr<GridType>(gridFactory.createGrid());

  // Reorder boundary IDs according to the inserction index
  const auto leafGridView(grid->leafGridView());
  if(!boundaryIDs.empty())
  {
    std::vector<int> tempIDs(boundaryIDs.size(),0);
    for(const auto& entity:elements(leafGridView))
      for(const auto& intersection:intersections(leafGridView,entity))
        if(intersection.boundary())
          tempIDs[intersection.boundarySegmentIndex()]=boundaryIDs[gridFactory.insertionIndex(intersection)];
    boundaryIDs=std::move(tempIDs);
  }

  // Load balancing and refinement
  grid->loadBalance();
  if ( refinements > 0 )
    grid->globalRefine( refinements );

  // Do some tests to make sure the grid has been properly read
  gridcheck(*grid);

  // Write MSH
  Dune::GmshWriter<typename GridType::LeafGridView> writer( leafGridView );
  writer.setPrecision(10);
  const std::string outputName("./"+gridName+"-"+gridManagerName+"-gmshtest-write.msh");
  writer.write(outputName);
  if(!elementsIDs.empty())
  {
    const std::string outputNameEntity("./"+gridName+"-"+gridManagerName+"-gmshtest-write-entity.msh");
    writer.write(outputNameEntity,elementsIDs);
  }
  if((!boundaryIDs.empty())&&(!elementsIDs.empty()))
  {
    const std::string outputNameBoundary("./"+gridName+"-"+gridManagerName+"-gmshtest-write-boundary.msh");
    writer.write(outputNameBoundary,elementsIDs,boundaryIDs);
  }

  // Write VTK
  std::ostringstream vtkName;
  vtkName << "./" << gridName << "-gmshtest-" << refinements;
  VTKWriter<typename GridType::LeafGridView> vtkWriter( leafGridView );
  vtkWriter.write( vtkName.str() );
  std::cout<<std::endl;
}


int main( int argc, char** argv )
try
{
  MPIHelper::instance( argc, argv );
  const int refinements = ( argc > 1 ) ? atoi( argv[1] ) : 0;
  const std::string path(static_cast<std::string>(DUNE_GRID_EXAMPLE_GRIDS_PATH)+"gmsh/");

#if GMSH_UGGRID
  testReadingAndWritingGrid<UGGrid<2> >( path, "curved2d", "UGGrid-2D", refinements );
  testReadingAndWritingGrid<UGGrid<2> >( path, "circle2ndorder", "UGGrid-2D", refinements );
  testReadingAndWritingGrid<UGGrid<2> >( path, "unitsquare_quads_2x2", "UGGrid-2D", refinements );
  testReadingAndWritingGrid<UGGrid<2> >( path, "hybrid-testgrid-2d", "UGGrid-2D", refinements );
  testReadingAndWritingGrid<UGGrid<3> >( path, "pyramid", "UGGrid-3D", refinements );
  testReadingAndWritingGrid<UGGrid<3> >( path, "pyramid2ndorder", "UGGrid-3D", refinements );
  testReadingAndWritingGrid<UGGrid<3> >( path, "hybrid-testgrid-3d", "UGGrid-3D", refinements );
#endif

#if GMSH_ALBERTAGRID
#if ALBERTA_DIM==2
  testReadingAndWritingGrid<AlbertaGrid<2> >( path, "curved2d", "AlbertaGrid-2D", refinements );
#endif
#if ALBERTA_DIM==3
  testReadingAndWritingGrid<AlbertaGrid<2> >( path, "sphere", "AlbertaGrid-2D", refinements );
  testReadingAndWritingGrid<AlbertaGrid<3> >( path, "pyramid", "AlbertaGrid-3D", refinements );
#endif
#endif

#if GMSH_ONEDGRID
  testReadingAndWritingGrid<OneDGrid>( path, "oned-testgrid", "OneDGrid", refinements );
#endif

  return 0;

}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e) {
  std::cerr << e.what() << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
