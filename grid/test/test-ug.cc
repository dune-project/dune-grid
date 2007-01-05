// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>

/*

   Instantiate UG-Grid and feed it to the generic gridcheck()

 */

#include <dune/grid/uggrid.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"

// Test parallel interface if a parallel UG is used
#ifdef ModelP
#include <mpi.h>
#include "checkcommunicate.cc"
#endif

class ArcOfCircle : public Dune::BoundarySegment<2>
{
public:

  ArcOfCircle(const Dune::FieldVector<double,2>& center, double radius,
              double fromAngle, double toAngle)
    : center_(center), radius_(radius), fromAngle_(fromAngle), toAngle_(toAngle)
  {}

  Dune::FieldVector<double,2> operator()(const Dune::FieldVector<double,1>& local) const {

    double angle = fromAngle_ + local[0]*(toAngle_ - fromAngle_);

    Dune::FieldVector<double,2> result = center_;
    result[0] += radius_ * std::cos(angle);
    result[1] += radius_ * std::sin(angle);

    return result;
  }

  Dune::FieldVector<double,2> center_;

  double radius_;

  double fromAngle_;

  double toAngle_;
};


void makeHalfCircleQuad(Dune::UGGrid<2>& grid, bool parametrization)
{
  using namespace Dune;

  grid.createBegin();

  // /////////////////////////////
  //   Create boundary segments
  // /////////////////////////////
  if (parametrization) {

    FieldVector<double,2> center(0);
    center[1] = 15;

    std::vector<unsigned int> vertices(2);

    vertices[0] = 1;  vertices[1] = 2;
    grid.insertBoundarySegment(vertices, new ArcOfCircle(center, 15, M_PI, M_PI*4/3));

    vertices[0] = 2;  vertices[1] = 3;
    grid.insertBoundarySegment(vertices, new ArcOfCircle(center, 15, M_PI*4/3, M_PI*5/3));

    vertices[0] = 3;  vertices[1] = 0;
    grid.insertBoundarySegment(vertices, new ArcOfCircle(center, 15, M_PI*5/3, M_PI*2));

  }

  // ///////////////////////
  //   Insert vertices
  // ///////////////////////
  FieldVector<double,2> pos;

  pos[0] = 15;  pos[1] = 15;
  grid.insertVertex(pos);

  pos[0] = -15; pos[1] = 15;
  grid.insertVertex(pos);

  pos[0] = -7.5; pos[1] = 2.00962;
  grid.insertVertex(pos);

  pos[0] = 7.5; pos[1] = 2.00962;
  grid.insertVertex(pos);
  // /////////////////
  // Insert elements
  // /////////////////

  std::vector<unsigned int> cornerIDs(4);
  cornerIDs[0] = 0;
  cornerIDs[1] = 1;
  cornerIDs[2] = 3;
  cornerIDs[3] = 2;

  grid.insertElement(GeometryType(GeometryType::cube,2), cornerIDs);

  // //////////////////////////////////////
  //   Finish initialization
  // //////////////////////////////////////
  grid.createEnd();

}

void make2DTestGrid(Dune::UGGrid<2>& grid)
{
  // Start grid creation
  grid.createBegin();

  // The list of grid vertex positions
  int numVertices = 16;

  double vertices[16][2] = {{0, 0},
                            {0.5, 0},
                            {0.5, 0.5},
                            {0, 0.5},
                            {0.25, 0},
                            {0.5, 0.25},
                            {0.25, 0.5},
                            {0, 0.25},
                            {0.25, 0.25},
                            {1, 0},
                            {1, 0.5},
                            {0.75, 0.25},
                            {1, 1},
                            {0.5, 1},
                            {0, 1},
                            {0.25, 0.75}};

  // Create the grid vertices
  for (int i=0; i<numVertices; i++) {
    Dune::FieldVector<double,2> pos;
    pos[0] = vertices[i][0];
    pos[1] = vertices[i][1];
    grid.insertVertex(pos);
  }

  // Create the triangle elements
  int numTriangles = 2;
  unsigned int triangles[2][3] = {{9, 10, 11},
                                  {15, 13, 14}};

  for (int i=0; i<numTriangles; i++) {
    std::vector<unsigned int> cornerIDs(3);
    for (int j=0; j<3; j++)
      cornerIDs[j] = triangles[i][j];
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),cornerIDs);
  }

  // Create the quadrilateral elements
  int numQuadrilaterals = 9;
  unsigned int quadrilaterals[9][4] = {{0, 4, 7, 8},
                                       {4, 1, 8, 5},
                                       {8, 5, 6, 2},
                                       {7, 8, 3, 6},
                                       {1, 9, 5, 11},
                                       {5, 11, 2, 10},
                                       {2, 10, 13, 12},
                                       {3, 6, 14, 15},
                                       {6, 2, 15, 13}};

  for (int i=0; i<numQuadrilaterals; i++) {
    std::vector<unsigned int> cornerIDs(4);
    for (int j=0; j<4; j++)
      cornerIDs[j] = quadrilaterals[i][j];
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2),cornerIDs);
  }

  // Finish initialization
  grid.createEnd();
}


void make3DTestGrid(Dune::UGGrid<3>& grid)
{
  // Start grid creation
  grid.createBegin();

  // The list of grid vertex positions
  int numVertices = 61;

  double vertices[61][3] = {{0, 0, 0},
                            {0.5, 0, 0},
                            {0.5, 0.5, 0},
                            {0, 0.5, 0},
                            {0, 0, 0.5},
                            {0.5, 0, 0.5},
                            {0.5, 0.5, 0.5},
                            {0, 0.5, 0.5},
                            {0.25, 0, 0},
                            {0.5, 0.25, 0},
                            {0.25, 0.5, 0},
                            {0, 0.25, 0},
                            {0, 0, 0.25},
                            {0.5, 0, 0.25},
                            {0.5, 0.5, 0.25},
                            {0, 0.5, 0.25},
                            {0.25, 0, 0.5},
                            {0.5, 0.25, 0.5},
                            {0.25, 0.5, 0.5},
                            {0, 0.25, 0.5},
                            {0.25, 0.25, 0},
                            {0.25, 0, 0.25},
                            {0.5, 0.25, 0.25},
                            {0.25, 0.5, 0.25},
                            {0, 0.25, 0.25},
                            {0.25, 0.25, 0.5},
                            {0.25, 0.25, 0.25},
                            {1, 0, 0},
                            {1, 0.5, 0},
                            {1, 0, 0.5},
                            {1, 0.5, 0.5},
                            {0.75, 0.25, 0.25},
                            {1, 1, 0},
                            {0.5, 1, 0},
                            {1, 1, 0.5},
                            {0.5, 1, 0.5},
                            {0.75, 0.75, 0.25},
                            {0, 1, 0},
                            {0, 1, 0.5},
                            {0.25, 0.75, 0.25},
                            {0, 0, 1},
                            {0.5, 0, 1},
                            {0.5, 0.5, 1},
                            {0, 0.5, 1},
                            {0.25, 0.25, 0.75},
                            {1, 0, 1},
                            {1, 0.5, 1},
                            {0.75, 0.25, 0.75},
                            {1, 1, 1},
                            {0.5, 1, 1},
                            {0, 1, 1},
                            {0.25, 0.75, 0.75},
                            {1.5, 0, 0},
                            {1.5, 0.5, 0},
                            {1.5, 1, 0},
                            {1.5, 0, 0.5},
                            {1.5, 0.5, 0.5},
                            {1.5, 1, 0.5},
                            {1.5, 0, 1},
                            {1.5, 0.5, 1},
                            {1.5, 1, 1}};

  // Create the grid vertices
  for (int i=0; i<numVertices; i++) {
    Dune::FieldVector<double,3> pos;
    for (int j=0; j<3; j++)
      pos[j] = vertices[i][j];
    grid.insertVertex(pos);
  }



  // Create the tetrahedron elements
  int numTetrahedra = 54;
  unsigned int tetrahedra[54][4] = {{10, 29, 3, 32},
                                    {10, 2, 28, 32},
                                    {10, 28, 29, 32},
                                    {14, 28, 2, 32},
                                    {14, 6, 30, 32},
                                    {14, 30, 28, 32},
                                    {15, 31, 7, 32},
                                    {15, 3, 29, 32},
                                    {15, 29, 31, 32},
                                    {18, 30, 6, 32},
                                    {18, 7, 31, 32},
                                    {18, 31, 30, 32},
                                    {15, 29, 3, 37},
                                    {15, 7, 31, 37},
                                    {15, 31, 29, 37},
                                    {15, 36, 7, 37},
                                    {15, 3, 34, 37},
                                    {15, 34, 36, 37},
                                    {11, 38, 4, 40},
                                    {11, 3, 34, 40},
                                    {11, 34, 38, 40},
                                    {15, 34, 3, 40},
                                    {15, 7, 36, 40},
                                    {15, 36, 34, 40},
                                    {16, 39, 8, 40},
                                    {16, 4, 38, 40},
                                    {16, 38, 39, 40},
                                    {19, 36, 7, 40},
                                    {19, 8, 39, 40},
                                    {19, 39, 36, 40},
                                    {17, 42, 6, 45},
                                    {17, 5, 41, 45},
                                    {17, 41, 42, 45},
                                    {18, 43, 7, 45},
                                    {18, 6, 42, 45},
                                    {18, 42, 43, 45},
                                    {19, 44, 8, 45},
                                    {19, 7, 43, 45},
                                    {19, 43, 44, 45},
                                    {20, 41, 5, 45},
                                    {20, 8, 44, 45},
                                    {20, 44, 41, 45},
                                    {18, 31, 7, 48},
                                    {18, 6, 30, 48},
                                    {18, 30, 31, 48},
                                    {18, 42, 6, 48},
                                    {18, 7, 43, 48},
                                    {18, 43, 42, 48},
                                    {19, 39, 8, 52},
                                    {19, 7, 36, 52},
                                    {19, 36, 39, 52},
                                    {19, 43, 7, 52},
                                    {19, 8, 44, 52},
                                    {19, 44, 43, 52}};

  for (int i=0; i<numTetrahedra; i++) {
    std::vector<unsigned int> cornerIDs(4);
    for (int j=0; j<4; j++)
      cornerIDs[j] = tetrahedra[i][j]-1;
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3),cornerIDs);
  }

  // Create the pyramid elements
  int numPyramids = 27;
  unsigned int pyramids[27][5] = {{28, 30, 31, 29, 32},
                                  {10, 23, 14, 2, 32},
                                  {14, 23, 18, 6, 32},
                                  {18, 23, 15, 7, 32},
                                  {15, 23, 10, 3, 32},
                                  {3, 29, 33, 34, 37},
                                  {29, 31, 35, 33, 37},
                                  {33, 35, 36, 34, 37},
                                  {7, 36, 35, 31, 37},
                                  {11, 24, 15, 3, 40},
                                  {15, 24, 19, 7, 40},
                                  {19, 24, 16, 8, 40},
                                  {16, 24, 11, 4, 40},
                                  {34, 36, 39, 38, 40},
                                  {20, 26, 19, 8, 45},
                                  {19, 26, 18, 7, 45},
                                  {18, 26, 17, 6, 45},
                                  {17, 26, 20, 5, 45},
                                  {41, 44, 43, 42, 45},
                                  {6, 42, 46, 30, 48},
                                  {30, 46, 47, 31, 48},
                                  {31, 47, 43, 7, 48},
                                  {42, 43, 47, 46, 48},
                                  {7, 43, 50, 36, 52},
                                  {36, 50, 51, 39, 52},
                                  {39, 51, 44, 8, 52},
                                  {44, 51, 50, 43, 52}};

  for (int i=0; i<numPyramids; i++) {
    std::vector<unsigned int> cornerIDs(5);
    for (int j=0; j<5; j++)
      cornerIDs[j] = pyramids[i][j]-1;
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,3),cornerIDs);
  }

  // Create the prism elements
  int numPrisms = 8;
  unsigned int prisms[8][6] = {{28, 53, 29, 30, 56, 31},
                               {53, 54, 29, 56, 57, 31},
                               {30, 56, 31, 46, 59, 47},
                               {56, 57, 31, 59, 60, 47},
                               {29, 54, 33, 31, 57, 35},
                               {54, 55, 33, 57, 58, 35},
                               {31, 57, 35, 47, 60, 49},
                               {57, 58, 35, 60, 61, 49}};


  for (int i=0; i<numPrisms; i++) {
    std::vector<unsigned int> cornerIDs(6);
    for (int j=0; j<6; j++)
      cornerIDs[j] = prisms[i][j]-1;
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::prism,3),cornerIDs);
  }

  // Create the hexahedron elements
  int numHexahedra = 9;
  unsigned int hexahedra[9][8] = {{1, 9, 12, 21, 13, 22, 25, 27},
                                  {9, 2, 21, 10, 22, 14, 27, 23},
                                  {21, 10, 11, 3, 27, 23, 24, 15},
                                  {12, 21, 4, 11, 25, 27, 16, 24},
                                  {13, 22, 25, 27, 5, 17, 20, 26},
                                  {22, 14, 27, 23, 17, 6, 26, 18},
                                  {27, 23, 24, 15, 26, 18, 19, 7},
                                  {25, 27, 16, 24, 20, 26, 8, 19},
                                  {7, 31, 36, 35, 43, 47, 50, 49}};


  for (int i=0; i<numHexahedra; i++) {
    std::vector<unsigned int> cornerIDs(8);
    for (int j=0; j<8; j++)
      cornerIDs[j] = hexahedra[i][j]-1;
    grid.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3),cornerIDs);
  }

  // Finish initialization
  grid.createEnd();
}


template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  typedef typename GridType::template Codim<0>::LeafIterator LeafIterator;

  int count = 0;

  LeafIterator endit = grid.template leafend  <0> ();
  for(LeafIterator it = grid.template leafbegin<0> (); it != endit ; ++it )
  {
    if(num == count) grid.mark( ref, it );
    count++;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}

int main (int argc , char **argv) try
{

#ifdef ModelP
  // initialize MPI
  MPI_Init(&argc,&argv);
#endif

  // ////////////////////////////////////////////////////////////////////////
  //  Do the standard grid test for a 2d UGGrid
  // ////////////////////////////////////////////////////////////////////////
  // extra-environment to check destruction
  {

    std::cout << "Testing UGGrid<2> and UGGrid<3> simultaneously" << std::endl;

    Dune::UGGrid<2> grid2d;
    Dune::UGGrid<3> grid3d;

    make2DTestGrid(grid2d);
    make3DTestGrid(grid3d);

    // check macro grid
    gridcheck(grid2d);
    gridcheck(grid3d);

    // create hybrid grid
    markOne(grid2d,0,1) ;
    markOne(grid3d,0,1) ;

    gridcheck(grid2d);
    gridcheck(grid3d);

    grid2d.globalRefine(1);
    grid3d.globalRefine(1);
    gridcheck(grid2d);
    gridcheck(grid3d);

    // check the method geometryInFather()
    checkGeometryInFather(grid2d);
    checkGeometryInFather(grid3d);

    // check the intersection iterator
    checkIntersectionIterator(grid2d);
    checkIntersectionIterator(grid3d);

#ifdef ModelP
    // check communication interface
    checkCommunication(grid2d,-1,Dune::dvverb);
    for(int l=0; l<=grid2d.maxLevel(); ++l)
      checkCommunication(grid2d,l,Dune::dvverb);
#endif
  }

  // ////////////////////////////////////////////////////////////////////////
  //   Check whether geometryInFather returns equal results with and
  //   without parametrized boundaries
  // ////////////////////////////////////////////////////////////////////////

  Dune::UGGrid<2> gridWithParametrization, gridWithoutParametrization;

  // make grids
  makeHalfCircleQuad(gridWithoutParametrization, false);
  makeHalfCircleQuad(gridWithParametrization, true);

  // make grids again just to check this is possible
#if 0
  makeHalfCircleQuad(gridWithoutParametrization, false);
  makeHalfCircleQuad(gridWithParametrization, true);
#endif

  gridWithParametrization.globalRefine(1);
  gridWithoutParametrization.globalRefine(1);

  typedef Dune::UGGrid<2>::Codim<0>::LevelIterator ElementIterator;
  ElementIterator eIt    = gridWithParametrization.lbegin<0>(1);
  ElementIterator eWoIt  = gridWithoutParametrization.lbegin<0>(1);
  ElementIterator eEndIt = gridWithParametrization.lend<0>(1);

  for (; eIt!=eEndIt; ++eIt, ++eWoIt) {

    // The grids where constructed identically and they are traversed identically
    // Thus their respective output from geometryInFather should be the same
    for (int i=0; i<eIt->geometry().corners(); i++) {

      Dune::FieldVector<double,2> diff = eIt->geometryInFather()[i] - eWoIt->geometryInFather()[i];

      if ( diff.two_norm() > 1e-5 )
        DUNE_THROW(Dune::GridError, "output of geometryInFather() depends on boundary parametrization!");

    }

  }

  // ////////////////////////////////////////////////////////////////////////
  //   Check whether copies of elements have the same global ID
  // ////////////////////////////////////////////////////////////////////////

  {
    std::cout << "Testing if copies of elements have the same globalID." << std::endl;
    Dune::UGGrid<2> locallyRefinedGrid;

    locallyRefinedGrid.setRefinementType(Dune::UGGrid<2>::COPY);

    typedef Dune::UGGrid<2>::Codim<0>::LevelIterator ElementIterator;
    typedef Dune::UGGrid<2>::Codim<0>::HierarchicIterator HierarchicIterator;
    typedef Dune::UGGrid<2>::Traits::GlobalIdSet GlobalIdSet;

    // make grids
    makeHalfCircleQuad(locallyRefinedGrid, false);

    markOne(locallyRefinedGrid,0,1);
    markOne(locallyRefinedGrid,0,1);

    const GlobalIdSet& globalIdSet = locallyRefinedGrid.globalIdSet();

    for (int level=0; level<locallyRefinedGrid.maxLevel(); ++level)
    {
      ElementIterator eIt = locallyRefinedGrid.lbegin<0>(level);
      ElementIterator eEnd = locallyRefinedGrid.lend<0>(level);
      for(; eIt!=eEnd; ++eIt)
      {
        int children = 0;
        GlobalIdSet::IdType globalChildId;

        HierarchicIterator hIt = eIt->hbegin(level+1);
        HierarchicIterator hEnd = eIt->hend(level+1);


        for( ; hIt!=hEnd; ++hIt)
        {
          globalChildId = globalIdSet.id<0>(*hIt);
          ++children;
        }
        if ((children == 1) && (globalIdSet.id<0>(*eIt) != globalChildId))
          DUNE_THROW(Dune::GridError, "Copy of element has different globalId!");
      }
    }
  }


#ifdef ModelP
  // Terminate MPI
  MPI_Finalize();
#endif

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
