// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/type_traits.hh>

// dune grid includes
#if HAVE_DUNE_UGGRID
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

#include <dune/common/exceptions.hh>
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

template<class T>
T &discarded(T &&v) { return v; }

template<class Grid, class... Args>
using read_gf_result_t =
  decltype(GmshReader<Grid>::read(std::declval<GridFactory<Grid>&>(),
                                  std::string{}, std::declval<Args>()...));

template <typename GridType>
void testReadingAndWritingGrid( const std::string& path, const std::string& gridName,
                                const std::string& gridManagerName, int refinements,
                                bool expectsBoundarySegments = false)
{
  // Read the grid
  std::cout<<"Using "<<gridManagerName<<std::endl;
  GridFactory<GridType> gridFactory;
  const std::string inputName(path+gridName+".msh");
  std::cout<<"Reading mesh file "<<inputName<<std::endl;
  auto reader = GmshReader<GridType>(inputName, gridFactory);
  auto grid = gridFactory.createGrid();
  auto elementsIDs = reader.extractElementData();
  auto boundaryIDs = reader.extractBoundaryData();

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
  else
    if (expectsBoundarySegments)
        DUNE_THROW(Dune::IOError, "Expected boundary segment markers found none!");

  // Load balancing and refinement
  grid->loadBalance();
  if ( refinements > 0 )
    grid->globalRefine( refinements );

  // Do some tests to make sure the grid has been properly read
  gridcheck(*grid);

  // Write MSH
  Dune::GmshWriter<typename GridType::LeafGridView> writer( leafGridView );
  writer.setPrecision(10);
  const std::string outputName(gridName+"-"+gridManagerName+"-gmshtest-write.msh");
  writer.write(outputName);
  if(!elementsIDs.empty())
  {
    const std::string outputNameEntity(gridName+"-"+gridManagerName+"-gmshtest-write-entity.msh");
    writer.write(outputNameEntity,elementsIDs);
  }
  if((!boundaryIDs.empty())&&(!elementsIDs.empty()))
  {
    const std::string outputNameBoundary(gridName+"-"+gridManagerName+"-gmshtest-write-boundary.msh");
    writer.write(outputNameBoundary,elementsIDs,boundaryIDs);
  }

  // Write VTK
  std::ostringstream vtkName;
  vtkName << gridName << "-gmshtest-" << refinements;
  VTKWriter<typename GridType::LeafGridView> vtkWriter( leafGridView );
  vtkWriter.write( vtkName.str() );
  std::cout<<std::endl;

  //
  // test all signatures of the read method
  //

  // Test whether grid can be read without giving the gridfactory explicitly
  {
    // unique_ptr
    std::unique_ptr<GridType> gridUnique =
      GmshReader<GridType>::read(inputName);
  }

  {
    // shared_ptr
    std::shared_ptr<GridType> gridShared =
      GmshReader<GridType>::read(inputName);
  }

  // test deprecated reading without gridfactory but with data
  {
    auto read = [&] (auto... args)
    {
      std::vector<int> boundaryData;
      std::vector<int> elementData;
      DUNE_NO_DEPRECATED_BEGIN
      auto gridp = GmshReader<GridType>::read(inputName, elementData,
                                              boundaryData, args...);
      DUNE_NO_DEPRECATED_END
      static_assert(std::is_same<std::remove_reference_t<decltype(*gridp)>,
                                 GridType>::value,
                    "GmshReader::read() return type is wrong");
    };

    read();
    read(false);
    read(false, false);
    read(false, true);
    read(true);
    read(true, false);
    read(true, true);
  }

  // test reading with gridfactory but without data
  {
    auto read = [&] (auto... args)
    {
      GridFactory<GridType> factory;
      GmshReader<GridType>::read(factory, inputName, args...);
    };

    read();
    read(false);
    read(false, false);
    read(false, true);
    read(true);
    read(true, false);
    read(true, true);
  }

  // test reading with gridfactory and data
  {
    auto read = [&] (auto &&... args)
    {
      return GmshReader<GridType>::read(discarded(GridFactory<GridType>{}),
                                        inputName,
                                        std::forward<decltype(args)>(args)...);
    };

    // boundary data provided, element data provided
    read(discarded(std::vector<int>{}), discarded(std::vector<int>{}));
    read(discarded(std::vector<int>{}), discarded(std::vector<int>{}), false);
    read(discarded(std::vector<int>{}), discarded(std::vector<int>{}), true);

    // boundary data provided, element data omitted
    read(discarded(std::vector<int>{}), std::ignore);
    read(discarded(std::vector<int>{}), std::ignore, false);
    read(discarded(std::vector<int>{}), std::ignore, true);

    // boundary data omitted, segment insertion omitted, element data provided
    read(std::ignore, discarded(std::vector<int>{}));
    read(std::ignore, discarded(std::vector<int>{}), false);
    read(std::ignore, discarded(std::vector<int>{}), true);

    // boundary data omitted, segment insertion omitted, element data omitted
    read(std::ignore, std::ignore);
    read(std::ignore, std::ignore, false);
    read(std::ignore, std::ignore, true);

    // boundary data omitted, segment insertion omitted, element data provided
    read(false, discarded(std::vector<int>{}));
    read(false, discarded(std::vector<int>{}), false);
    read(false, discarded(std::vector<int>{}), true);

    // boundary data omitted, segment insertion omitted, element data omitted
    read(false, std::ignore);
    read(false, std::ignore, false);
    read(false, std::ignore, true);

    // boundary data omitted, segment insertion done, element data provided
    read(true, discarded(std::vector<int>{}));
    read(true, discarded(std::vector<int>{}), false);
    read(true, discarded(std::vector<int>{}), true);

    // boundary data omitted, segment insertion done, element data omitted
    read(true, std::ignore);
    read(true, std::ignore, false);
    read(true, std::ignore, true);

    // check that certain signatures are rejected
    static_assert(!Std::is_detected_v<read_gf_result_t, GridType,
                                      std::vector<int>, std::vector<int>&>,
                  "rvalue boundary data should be rejected");
    static_assert(!Std::is_detected_v<read_gf_result_t, GridType,
                                      std::vector<int>&, std::vector<int>>,
                  "rvalue element data should be rejected");
  }

  // test reading with gridfactory and data (deprecated signatures)
  {
    auto read = [&] (bool verbose, bool insertBoundarySegments)
    {
      std::vector<int> boundaryData;
      std::vector<int> elementData;
      GridFactory<GridType> factory;
      DUNE_NO_DEPRECATED_BEGIN
      GmshReader<GridType>::read(factory, inputName, boundaryData, elementData,
                                 verbose, insertBoundarySegments);
      DUNE_NO_DEPRECATED_END
    };

    read(false, false);
    read(false, true);
    read(true, false);
    read(true, true);
  }
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
  testReadingAndWritingGrid<UGGrid<3> >( path, "unitcube", "UGGrid-3D", refinements );
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
  testReadingAndWritingGrid<OneDGrid>( path, "oned-testgrid", "OneDGrid", refinements, true);
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
