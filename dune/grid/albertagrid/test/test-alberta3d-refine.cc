// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif

// single unit tet
template<typename Grid>
class UnitTetMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::shared_ptr<Grid> create() {
    Dune::GridFactory<Grid> gf;

    // insert vertices
    Dune::FieldVector<typename Grid::ctype, 3> pos;
    pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
    pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);


    // insert elements
    std::vector<unsigned int> vid(4);
    vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);

    return std::shared_ptr<Grid>(gf.createGrid());
  }
};

// kuhn triangulation with 6 tets
template<typename Grid>
class KuhnTriangulatedUnitCubeMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::shared_ptr<Grid> create() {
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

    std::vector<unsigned int> vid(4);
    vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);

    gf.markLongestEdge();
    return std::shared_ptr<Grid>(gf.createGrid());
  }
};

// minimal triangulation with 5 tets, contains unit tet
template<typename Grid>
class MinTriangulatedUnitCubeMaker {
  static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
  static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
  static std::shared_ptr<Grid> create() {
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

    std::vector<unsigned int> vid(4);
    // tet at vertex 0
    vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 4; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    // tet at vertex 3
    vid[0] = 1; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    // central tet
    vid[0] = 1; vid[1] = 2; vid[2] = 4; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    // tet at vertex 5
    vid[0] = 1; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);
    // tet at vertex 6
    vid[0] = 2; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(Dune::GeometryTypes::tetrahedron, vid);

    gf.markLongestEdge();
    return std::shared_ptr<Grid>(gf.createGrid());
  }
};


int main(int /* argc */, char** /* argv*/)
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
      typedef Dune::AlbertaGrid<3, 3> Grid;

      std::cout << "The recursive-bisection refinement algorithm of alberta "
                << "cannot be used on arbitrary meshes (see Flyspry#569 for a "
                << "more comprehensive explanation).  However, the heuristic "
                << "used in the GridFactory can be tuned to support certain "
                << "types of meshes.  This test makes sure that meshes that "
                << "worked at one time in the past continue to work in the "
                << "future.  If a certain mesh does not work, this program "
                << "will usually generate a segmentation fault.\n"
                << std::endl;

      std::cout << "Checking unit tetrahedron..." << std::endl;
      std::cout << "Note: The unit tetrahdron check is a safety measure.  If "
                << "this test already produces a segfault, the problem is "
                << "probably not Albertas recursive bisection algorithm but "
                << "something else." << std::endl;
      UnitTetMaker<Grid>::create()->globalRefine(2);
      std::cout << "Checking unit tetrahedron: success\n" << std::endl;

      std::cout << "Checking 6-triangulation of unit cube..." << std::endl;
      KuhnTriangulatedUnitCubeMaker<Grid>::create()->globalRefine(2);
      std::cout << "Checking 6-triangulation of unit cube: success\n"
                << std::endl;

      std::cout << "Checking 5-triangulation of unit cube..." << std::endl;
      MinTriangulatedUnitCubeMaker<Grid>::create()->globalRefine(2);
      std::cout << "Checking 5-triangulation of unit cube: success\n"
                << std::endl;

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
