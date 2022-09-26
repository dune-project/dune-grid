// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief A unit test for the PersistentContainer
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;

template <int dw>
struct Data
{
  typedef FieldVector<double,dw> Coord;
  Data() : used(false) {}
  Data(const Coord &x) : coord(x), used(true) {}
  Coord coord;
  bool used;
};

template <class GridType>
bool test(GridType &grid)
{
  bool ret = true;
  typedef Data<GridType::dimensionworld> DataType;
  PersistentContainer<GridType,DataType> container0(grid,0);
  PersistentContainer<GridType,DataType> container1(grid,1);
  PersistentContainer<GridType,DataType> container2(grid,2);
  PersistentContainer<const GridType,DataType> ccontainer0(grid,0);
  PersistentContainer<const GridType,DataType> ccontainer1(grid,1);
  PersistentContainer<const GridType,DataType> ccontainer2(grid,2);

  typedef typename GridType::LeafGridView GridView;
  const GridView &view = grid.leafGridView();
  typedef typename GridView::template Codim<0>::Iterator EIterator;

  {
    const EIterator &eend = view.template end<0>();
    for(EIterator eit = view.template begin<0>(); eit != eend; ++eit)
    {
      auto geometry = eit->geometry();
      ccontainer0[*eit] = container0[*eit] = geometry.center();
      auto refElement = referenceElement(geometry);
      for (unsigned int i=0; i<eit->subEntities(1); ++i)
        ccontainer1(*eit,i) = container1(*eit,i) = geometry.global( refElement.position(i,1) );
      for (unsigned int i=0; i<eit->subEntities(2); ++i)
        container2(*eit,i) = geometry.global( refElement.position(i,2) );
    }
  }

  (void) container0.size();
  (void) container1.size();
  (void) container2.size();

  (void) ccontainer0.size();
  (void) ccontainer1.size();
  (void) ccontainer2.size();

  grid.globalRefine(1);
  container0.resize();
  container1.resize();
  container2.resize();
  ccontainer0.resize();
  ccontainer1.resize();
  ccontainer2.resize();

  {
    const EIterator &eend = view.template end<0>();
    for(EIterator eit = view.template begin<0>(); eit != eend; ++eit)
    {
      if (container0[*eit].used == true)
      {
        std::cout << "ERROR: a new element is marked as 'used' in the container - stop testing" << std::endl;
        ret = false;
        break;
      }
      if (ccontainer0[*eit].used == true)
      {
        std::cout << "ERROR: a new element is marked as 'used' in the const container - stop testing" << std::endl;
        ret = false;
        break;
      }
      typename GridType::template Codim<0>::Entity up = eit->father();
      while ( !container0[up].used )
      {
        if (up.level() == 0)
        {
          std::cout << "ERROR: could not find a father element in container - stop testing" << std::endl;
          ret = false;
          break;
        }
        up = up.father();
      }
      if ( ( container0[up].coord - up.geometry().center() ).two_norm() > 1e-8 )
      {
        std::cout << "ERROR: wrong data stored in container0 - stop testing" << std::endl;
        ret = false;
        break;
      }
      auto refElement = Dune::referenceElement< typename GridType::ctype, GridType::dimension >( eit->type() );
      for (unsigned int i=0; i<eit->subEntities(1); ++i)
        if ( ( container1(up,i).coord - up.geometry().global( refElement.position(i,1) ) ).two_norm() > 1e-8 )
        {
          std::cout << "ERROR: wrong data stored in container1 - stop testing" << std::endl;
          ret = false;
          break;
        }
      for (unsigned int i=0; i<eit->subEntities(2); ++i)
        if ( ( container2(up,i).coord - up.geometry().global( refElement.position(i,2) ) ).two_norm() > 1e-8 )
        {
          std::cout << "ERROR: wrong data stored in container1 - stop testing" << std::endl;
          ret = false;
          break;
        }
    }
  }
  return ret;
}

int main (int argc , char **argv)
try {

  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  // /////////////////////////////////////////////////////////////////////////////
  //   Test YaspGrid
  // /////////////////////////////////////////////////////////////////////////////
  {
    typedef YaspGrid<2> GridType;
    FieldVector<double,2> Len; Len = 1.0;
    std::array<int,2> s = { {2, 6} };
    std::bitset<2> p;
    int overlap = 1;
    GridType grid(Len,s,p,overlap);
    std::cout << "Testing YaspGrid" << std::endl;
    test(grid);
  }

#if HAVE_DUNE_UGGRID
  // /////////////////////////////////////////////////////////////////////////////
  //   Test UGGrid
  // /////////////////////////////////////////////////////////////////////////////
  {
    typedef UGGrid<2> GridType;
    FieldVector<double,2> lowerLeft{0.0, 0.0};
    FieldVector<double,2> upperRight{1.0, 1.0};
    std::array<unsigned int, 2> cells = { {2, 6} };
    const auto gridPtr = StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, cells);
    std::cout << "Testing UGGrid" << std::endl;
    test(*gridPtr);
  }
#endif

  return 0;
}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
