// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>

#include <dune/grid/sgrid.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkpartition.cc"

template<int d, int w>
void runtest()
{
  int n[] = { 5, 5, 5, 5 };
  double h[] = { 1.0, 2.0, 3.0, 4.0 };

  std::cout << std::endl << "SGrid<" << d << "," << w << ">" << std::endl;
  Dune::SGrid<d,w> g(n, h);
  gridcheck(g);

  g.globalRefine(1);
  checkGeometryInFather(g);
  checkIntersectionIterator(g);
  checkPartitionType( g.leafView() );

  std::cout << std::endl;
}

void testFS940()
{
  typedef Dune::SGrid<2,2> GridType;
  int n[2];
  double h[2];

  for (int i=0; i<2; ++i) {
    n[i] = 2;
    h[i] = 1.0;
  }

  GridType grid(n,h);
  typedef GridType::Codim<0>::LevelIterator ElementIterator;

  ElementIterator et = grid.lbegin<0>(0);
  ElementIterator et2 = grid.lbegin<0>(0);
  ++et2;

  typedef  GridType::LeafGridView GridView;
  typedef  GridView::IntersectionIterator IntersectionIterator;
  const GridView gv = grid.leafView();

  IntersectionIterator i1 = gv.ibegin(*et);
  IntersectionIterator i2 = gv.ibegin(*et2);

  // expected results:
  int res1[] = {0,1,1};
  int res2[] = {1,0,0};

  std::cout<<"i1: \n";
  std::cout<<" has neighbor "<<i1->neighbor()<<std::endl;
  std::cout<<" boundary "<<i1->boundary()<<std::endl;
  std::cout<<" boundaryId "<<i1->boundaryId()<<std::endl;

  if (i1->neighbor() != res1[0]
      || i1->boundary() != res1[1]
      || i1->boundaryId() != res1[2])
  {
    DUNE_THROW(Dune::Exception,
               "Wrong intersection information for i1: "
               << " has neighbor/boundary/boundaryId "
               << i1->neighbor() << "/" << i1->boundary() << "/" << i1->boundaryId()
               << "expected " << res1[0] << "/" << res1[1] << "/" << res1[2]);
  }

  std::cout<<"i2: \n";
  std::cout<<" has neighbor "<<i2->neighbor()<<std::endl;
  std::cout<<" boundary "<<i2->boundary()<<std::endl;
  std::cout<<" boundaryId "<<i2->boundaryId()<<std::endl;

  if (i2->neighbor() != res2[0]
      || i2->boundary() != res2[1]
      || i2->boundaryId() != res2[2])
  {
    DUNE_THROW(Dune::Exception,
               "Wrong intersection information for i2: "
               << " has neighbor/boundary/boundaryId "
               << i2->neighbor() << "/" << i2->boundary() << "/" << i2->boundaryId()
               << "expected " << res2[0] << "/" << res2[1] << "/" << res2[2]);
  }

  i1=i2;

  std::cout<<"i1 after i1=i2: \n";
  std::cout<<" has neighbor "<<i1->neighbor()<<std::endl;
  std::cout<<" boundary "<<i1->boundary()<<std::endl;
  std::cout<<" boundaryId "<<i1->boundaryId()<<std::endl;

  if (i1->neighbor() != res2[0]
      || i1->boundary() != res2[1]
      || i1->boundaryId() != res2[2])
  {
    DUNE_THROW(Dune::Exception,
               "Wrong intersection information for i1 after assignment: "
               << " has neighbor/boundary/boundaryId "
               << i1->neighbor() << "/" << i1->boundary() << "/" << i1->boundaryId()
               << "expected " << res2[0] << "/" << res2[1] << "/" << res2[2]);
  }

}

int main () {
  try {
    runtest<1,1>();
    runtest<2,2>();
    runtest<3,3>();
    //    runtest<4,4>();
    runtest<1,3>();
    runtest<2,3>();
    testFS940();
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
