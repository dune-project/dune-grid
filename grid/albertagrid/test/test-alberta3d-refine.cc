// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif

// minimal triangulation with 5 tets, contains unit tet
template<typename Grid>
class TriangulatedUnitCubeMaker {
  dune_static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  dune_static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static Dune::SmartPointer<Grid> create() {
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

    // tet at vertex 0
    vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 4; gf.insertElement(type, vid);
    // tet at vertex 3
    vid[0] = 1; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
    // central tet
    vid[0] = 1; vid[1] = 2; vid[2] = 4; vid[3] = 7; gf.insertElement(type, vid);
    // tet at vertex 5
    vid[0] = 1; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
    // tet at vertex 6
    vid[0] = 2; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);

    return gf.createGrid();
  }
};


int main(int argc, char** argv)
{

  try{
    // default exitcode 77 (=skipped); returned in case none of the supported
    // Grids were found
    int result = 77;

#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
    {
      std::cout << "The following stuff is expected to fail currently, because the recursive" << std::endl
                << "bisection algorithm of alberta will run into an endless recursion for" << std::endl
                << "certain grids.  See Flyspry#569." << std::endl;
      typedef Dune::AlbertaGrid<3, 3> Grid;
      Dune::SmartPointer<Grid> grid = TriangulatedUnitCubeMaker<Grid>::create();
      grid->globalRefine(2);

      result = 0;
    }
#endif

    return result;
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
