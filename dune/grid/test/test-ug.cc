// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>

/*

   Instantiate UG-Grid and feed it to the generic gridcheck()

 */

#include <dune/grid/uggrid.hh>
#include <doc/grids/gridfactory/hybridtestgrids.hh>

#include "gridcheck.cc"
#include "checkcommunicate.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"

#include <dune/common/parallel/mpihelper.hh>

using namespace Dune;

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


void makeHalfCircleQuad(Dune::UGGrid<2>& grid, bool boundarySegments, bool parametrization)
{
  Dune::GridFactory<Dune::UGGrid<2> > factory(&grid);

  // /////////////////////////////
  //   Create boundary segments
  // /////////////////////////////
  if (boundarySegments) {

    FieldVector<double,2> center(0);
    center[1] = 15;

    std::vector<unsigned int> vertices(2);

    if (parametrization) {

      vertices[0] = 1;  vertices[1] = 2;
      factory.insertBoundarySegment(vertices, shared_ptr<BoundarySegment<2,2> >(new ArcOfCircle(center, 15, M_PI, M_PI*4/3)));

      vertices[0] = 2;  vertices[1] = 3;
      factory.insertBoundarySegment(vertices, shared_ptr<BoundarySegment<2,2> >(new ArcOfCircle(center, 15, M_PI*4/3, M_PI*5/3)));

      vertices[0] = 3;  vertices[1] = 0;
      factory.insertBoundarySegment(vertices, shared_ptr<BoundarySegment<2,2> >(new ArcOfCircle(center, 15, M_PI*5/3, M_PI*2)));

    } else {

      vertices[0] = 1;  vertices[1] = 2;
      factory.insertBoundarySegment(vertices);

      vertices[0] = 2;  vertices[1] = 3;
      factory.insertBoundarySegment(vertices);

      vertices[0] = 3;  vertices[1] = 0;
      factory.insertBoundarySegment(vertices);

    }

  }

  // ///////////////////////
  //   Insert vertices
  // ///////////////////////
  FieldVector<double,2> pos;

  pos[0] = 15;  pos[1] = 15;
  factory.insertVertex(pos);

  pos[0] = -15; pos[1] = 15;
  factory.insertVertex(pos);

  pos[0] = -7.5; pos[1] = 2.00962;
  factory.insertVertex(pos);

  pos[0] = 7.5; pos[1] = 2.00962;
  factory.insertVertex(pos);

  // /////////////////
  // Insert elements
  // /////////////////

  std::vector<unsigned int> cornerIDs(4);
  cornerIDs[0] = 0;
  cornerIDs[1] = 1;
  cornerIDs[2] = 3;
  cornerIDs[3] = 2;

  factory.insertElement(GeometryType(GeometryType::cube,2), cornerIDs);

  // //////////////////////////////////////
  //   Finish initialization
  // //////////////////////////////////////
  factory.createGrid();

}


template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  typedef typename GridType::template Codim<0>::LeafIterator LeafIterator;

  int count = 0;

  LeafIterator endit = grid.template leafend  <0> ();
  for(LeafIterator it = grid.template leafbegin<0> (); it != endit ; ++it )
  {
    if(num == count) grid.mark( ref, *it );
    count++;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}

void generalTests(bool greenClosure)
{
  // /////////////////////////////////////////////////////////////////
  //   Prelude: the UGGrid implementation relies on the face that std::array,
  //   Dune::FieldVector, and C-arrays have the same memory layout.  Here
  //   we make sure this is really so.
  // /////////////////////////////////////////////////////////////////

  double cArray[3] = {1,2,3};

  for (int i=0; i<3; i++) {
    assert(cArray[i] == (*((FieldVector<double,3>*)&cArray))[i]);
    assert(cArray[i] == (*((array<double,3>*)&cArray))[i]);
  }

  // //////////////////////////////////////////////////////////
  //   Make some grids for testing
  // //////////////////////////////////////////////////////////

  std::auto_ptr<Dune::UGGrid<2> > grid2d(make2DHybridTestGrid<Dune::UGGrid<2> >());
  std::auto_ptr<Dune::UGGrid<3> > grid3d(make3DHybridTestGrid<Dune::UGGrid<3> >());

  // Switch of the green closure, if requested
  if (!greenClosure) {
    grid2d->setClosureType(UGGrid<2>::NONE);
    grid3d->setClosureType(UGGrid<3>::NONE);
  }

  // check macro grid
  gridcheck(*grid2d);
  gridcheck(*grid3d);

  // check communication interface
  checkCommunication(*grid2d,-1,Dune::dvverb);
  checkCommunication(*grid3d,-1,Dune::dvverb);
  for(int l=0; l<=grid2d->maxLevel(); ++l)
    checkCommunication(*grid2d,l,Dune::dvverb);
  for(int l=0; l<=grid3d->maxLevel(); ++l)
    checkCommunication(*grid3d,l,Dune::dvverb);

  // create hybrid grid
  markOne(*grid2d,0,1) ;
  markOne(*grid3d,0,1) ;

  gridcheck(*grid2d);
  gridcheck(*grid3d);

  // check communication interface
  checkCommunication(*grid2d,-1,Dune::dvverb);
  checkCommunication(*grid3d,-1,Dune::dvverb);
  for(int l=0; l<=grid2d->maxLevel(); ++l)
    checkCommunication(*grid2d,l,Dune::dvverb);
  for(int l=0; l<=grid3d->maxLevel(); ++l)
    checkCommunication(*grid3d,l,Dune::dvverb);

  grid2d->globalRefine(1);
  grid3d->globalRefine(1);

  gridcheck(*grid2d);
  gridcheck(*grid3d);
  // check communication interface
  checkCommunication(*grid2d,-1,Dune::dvverb);
  checkCommunication(*grid3d,-1,Dune::dvverb);
  for(int l=0; l<=grid2d->maxLevel(); ++l)
    checkCommunication(*grid2d,l,Dune::dvverb);
  for(int l=0; l<=grid3d->maxLevel(); ++l)
    checkCommunication(*grid3d,l,Dune::dvverb);

  // check geometry lifetime
  checkGeometryLifetime( grid2d->leafView() );
  checkGeometryLifetime( grid3d->leafView() );

  // check the method geometryInFather()
  checkGeometryInFather(*grid2d);
  checkGeometryInFather(*grid3d);

  // check the intersection iterator
  checkIntersectionIterator(*grid2d);
  checkIntersectionIterator(*grid3d);

}

int main (int argc , char **argv) try
{
  // use MPI helper to initialize MPI
  MPIHelper :: instance( argc, argv );

  // ////////////////////////////////////////////////////////////////////////
  //  Do the standard grid test for a 2d and a 3d UGGrid
  // ////////////////////////////////////////////////////////////////////////

  // Do the general tests for red/green refinement
  std::cout << "Testing UGGrid<2> and UGGrid<3> with red/green refinement" << std::endl;
  generalTests(true);

  // Do the general tests for nonconforming refinement
  std::cout << "Testing UGGrid<2> and UGGrid<3> with nonconforming refinement" << std::endl;
  generalTests(false);

  // ////////////////////////////////////////////////////////////////////////////
  //   Test whether I can create a grid with explict boundary segment ordering,
  //   but not parametrization functions (only 2d, so far)
  // ////////////////////////////////////////////////////////////////////////////
#ifdef ModelP
  {  // Force destruction of test grid right after the next test,
     // because for parallel UG there can be only one (2d-)grid at once.
#endif
  Dune::UGGrid<2> gridWithOrderedBoundarySegments;
  makeHalfCircleQuad(gridWithOrderedBoundarySegments, true, false);

  gridWithOrderedBoundarySegments.globalRefine(1);
  gridcheck(gridWithOrderedBoundarySegments);
#ifdef ModelP
  }
#endif
  // ////////////////////////////////////////////////////////////////////////
  //   Check whether geometryInFather returns equal results with and
  //   without parametrized boundaries
  // ////////////////////////////////////////////////////////////////////////

  // Only the sequential UG can provide more than one 2d- or 3d-grid at once.
  // Therefore we do not perform the following test for parallel UGGrid.
#if ! defined ModelP
  Dune::UGGrid<2> gridWithParametrization, gridWithoutParametrization;

  // make grids
  makeHalfCircleQuad(gridWithoutParametrization, false, false);
  makeHalfCircleQuad(gridWithParametrization, true, true);

  // make grids again just to check this is possible
  makeHalfCircleQuad(gridWithoutParametrization, false, false);
  makeHalfCircleQuad(gridWithParametrization, true, true);

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

      Dune::FieldVector<double,2> diff = eIt->geometryInFather().corner(i) - eWoIt->geometryInFather().corner(i);

      if ( diff.two_norm() > 1e-5 )
        DUNE_THROW(Dune::GridError, "output of geometryInFather() depends on boundary parametrization!");

    }

  }
#endif
  // ////////////////////////////////////////////////////////////////////////
  //   Check whether copies of elements have the same global ID
  // ////////////////////////////////////////////////////////////////////////

  {
    std::cout << "Testing if copies of elements have the same id." << std::endl;
    Dune::UGGrid<2> locallyRefinedGrid;

    locallyRefinedGrid.setRefinementType(Dune::UGGrid<2>::COPY);

    typedef Dune::UGGrid<2>::Codim<0>::LevelIterator ElementIterator;
    typedef Dune::UGGrid<2>::HierarchicIterator HierarchicIterator;
    typedef Dune::UGGrid<2>::Traits::GlobalIdSet GlobalIdSet;
    typedef Dune::UGGrid<2>::Traits::LocalIdSet LocalIdSet;

    // make grids
    makeHalfCircleQuad(locallyRefinedGrid, false, false);

    markOne(locallyRefinedGrid,0,1);
    markOne(locallyRefinedGrid,0,1);

    const GlobalIdSet& globalIdSet = locallyRefinedGrid.globalIdSet();
    const LocalIdSet&  localIdSet  = locallyRefinedGrid.localIdSet();

    for (int level=0; level<locallyRefinedGrid.maxLevel(); ++level)
    {
      ElementIterator eIt = locallyRefinedGrid.lbegin<0>(level);
      ElementIterator eEnd = locallyRefinedGrid.lend<0>(level);
      for(; eIt!=eEnd; ++eIt)
      {
        int children = 0;
        GlobalIdSet::IdType globalChildId;
        LocalIdSet::IdType localChildId;

        HierarchicIterator hIt = eIt->hbegin(level+1);
        HierarchicIterator hEnd = eIt->hend(level+1);


        for( ; hIt!=hEnd; ++hIt)
        {
          globalChildId = globalIdSet.id<0>(*hIt);
          localChildId =  localIdSet.id<0>(*hIt);
          ++children;
        }
        if (children != 1)
          continue;

        if (globalIdSet.id<0>(*eIt) != globalChildId)
          DUNE_THROW(Dune::GridError, "Copy of element has different globalId!");

        if (localIdSet.id<0>(*eIt) != localChildId)
          DUNE_THROW(Dune::GridError, "Copy of element has different localId!");
      }
    }
  }

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
