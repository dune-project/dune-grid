// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

template<typename Grid>
Grid* createTriangulatedCube() {
  Dune::GridFactory<Grid> gf;

  Dune::FieldVector<typename Grid::ctype, 3> pos;

  pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
  pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

  Dune::GeometryType type;
  type.makeTetrahedron();
  std::vector<unsigned int> vid(4);

  vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
  vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
  vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
  vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
  vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
  vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);

  return gf.createGrid();
}

template<typename Grid>
Grid* createTetrahedron() {
  Dune::GridFactory<Grid> gf;

  Dune::FieldVector<typename Grid::ctype, 3> pos;

  pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
  pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

  Dune::GeometryType type;
  type.makeTetrahedron();
  std::vector<unsigned int> vid(4);

  vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3; gf.insertElement(type, vid);

  return gf.createGrid();
}


//===============================================================
// Main program with grid setup
//===============================================================

typedef Dune::UGGrid<3> Grid;

int main(int argc, char** argv)
{
  // refine this many times before using the grid
  unsigned level = 1;

  // Calculate timestep
  try{
    Dune::shared_ptr<Grid> grid(createTriangulatedCube<Grid>());
    grid->globalRefine(level);

    // Dune::VTKWriter<Grid::LeafGridView> writer(grid->leafView());
    // writer.write("grid");

    return 0;
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
}
