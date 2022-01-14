// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>

/*
   Instantiate UG-Grid and feed it to the generic gridcheck()
 */
#include <dune/grid/uggrid.hh>
#include <doc/grids/gridfactory/hybridtestgrids.hh>

#include "gridcheck.hh"
#include "checkcommunicate.hh"
#include "checkgeometryinfather.hh"
#include "checkintersectionit.hh"
#include "checkpartition.hh"


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

    FieldVector<double,2> center = {0,15};

    if (parametrization) {

      factory.insertBoundarySegment({1,2}, std::make_shared<ArcOfCircle>(center, 15, M_PI, M_PI*4/3));
      factory.insertBoundarySegment({2,3}, std::make_shared<ArcOfCircle>(center, 15, M_PI*4/3, M_PI*5/3));
      factory.insertBoundarySegment({3,0}, std::make_shared<ArcOfCircle>(center, 15, M_PI*5/3, M_PI*2));

    } else {

      factory.insertBoundarySegment({1,2});
      factory.insertBoundarySegment({2,3});
      factory.insertBoundarySegment({3,0});

    }

  }

  // ///////////////////////
  //   Insert vertices
  // ///////////////////////

  factory.insertVertex({15,15});
  factory.insertVertex({-15,15});
  factory.insertVertex({-7.5,2.00962});
  factory.insertVertex({7.5,2.00962});

  // /////////////////
  // Insert elements
  // /////////////////

  factory.insertElement(GeometryTypes::quadrilateral, {0,1,3,2});

  // //////////////////////////////////////
  //   Finish initialization
  // //////////////////////////////////////

  // The factory returns the grid object we have given it as a std::unique_ptr.
  // Release it, to keep the std::unique_ptr destructor from killing it at
  // the end of this very method.
  // Yes, that is a bit weird.  But in particular the calling method tests whether
  // the GridFactory can reuse an existing UGGrid object, and I would like to keep
  // that for the time being.
  factory.createGrid().release();

}


template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  int count = 0;

  for(const auto& element : elements(grid.leafGridView()))
  {
    if(num == count) grid.mark( ref, element );
    count++;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}

void generalTests(bool greenClosure)
{
  // /////////////////////////////////////////////////////////////////
  //   Prelude: the UGGrid implementation relies on the fact that std::array,
  //   Dune::FieldVector, and C-arrays have the same memory layout.  Here
  //   we make sure this is really so.
  // /////////////////////////////////////////////////////////////////

#ifndef NDEBUG
  double cArray[3] = {1,2,3};

  for (int i=0; i<3; i++) {
    assert(cArray[i] == (*((FieldVector<double,3>*)&cArray))[i]);
    assert(cArray[i] == (*((std::array<double,3>*)&cArray))[i]);
  }
#endif

  // //////////////////////////////////////////////////////////
  //   Make some grids for testing
  // //////////////////////////////////////////////////////////

  std::unique_ptr<Dune::UGGrid<2> > grid2d(make2DHybridTestGrid<Dune::UGGrid<2> >());
  std::unique_ptr<Dune::UGGrid<3> > grid3d(make3DHybridTestGrid<Dune::UGGrid<3> >());

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
  checkGeometryLifetime( grid2d->leafGridView() );
  checkGeometryLifetime( grid3d->leafGridView() );

  // check the method geometryInFather()
  checkGeometryInFather(*grid2d);
  checkGeometryInFather(*grid3d);

  // check the intersection iterator
  checkIntersectionIterator(*grid2d);
  checkIntersectionIterator(*grid3d);

  // Check partition iterators
  checkPartitionType( grid2d->leafGridView() );
  for( int i = 0; i <= grid2d->maxLevel(); ++i )
    checkPartitionType( grid2d->levelGridView( i ) );
  checkPartitionType( grid3d->leafGridView() );
  for( int i = 0; i <= grid3d->maxLevel(); ++i )
    checkPartitionType( grid3d->levelGridView( i ) );

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

  Dune::UGGrid<2> gridWithOrderedBoundarySegments;
  makeHalfCircleQuad(gridWithOrderedBoundarySegments, true, false);

  gridWithOrderedBoundarySegments.globalRefine(1);
  gridcheck(gridWithOrderedBoundarySegments);

  // ////////////////////////////////////////////////////////////////////////
  //   Check whether geometryInFather returns equal results with and
  //   without parametrized boundaries
  // ////////////////////////////////////////////////////////////////////////

  // Only the sequential UG can provide more than one 2d- or 3d-grid at once.
  // Therefore we do not perform the following test for parallel UGGrid.

  Dune::UGGrid<2> gridWithParametrization, gridWithoutParametrization;

  // make grids
  makeHalfCircleQuad(gridWithoutParametrization, false, false);
  makeHalfCircleQuad(gridWithParametrization, true, true);

  // make grids again just to check this is possible
  makeHalfCircleQuad(gridWithoutParametrization, false, false);
  makeHalfCircleQuad(gridWithParametrization, true, true);

  gridWithParametrization.globalRefine(1);
  gridWithoutParametrization.globalRefine(1);

  auto eIt    = gridWithParametrization.levelGridView(1).begin<0>();
  auto eWoIt  = gridWithoutParametrization.levelGridView(1).begin<0>();
  auto eEndIt = gridWithParametrization.levelGridView(1).end<0>();

  for (; eIt!=eEndIt; ++eIt, ++eWoIt) {

    // The grids where constructed identically and they are traversed identically
    // Thus their respective output from geometryInFather should be the same
    for (int i=0; i<eIt->geometry().corners(); i++) {

      Dune::FieldVector<double,2> diff = eIt->geometryInFather().corner(i) - eWoIt->geometryInFather().corner(i);

      if ( diff.two_norm() > 1e-5 )
        DUNE_THROW(Dune::GridError, "output of geometryInFather() depends on boundary parametrization!");

    }

  }

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
      ElementIterator eIt = locallyRefinedGrid.levelGridView(level).begin<0>();
      ElementIterator eEnd = locallyRefinedGrid.levelGridView(level).end<0>();
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
