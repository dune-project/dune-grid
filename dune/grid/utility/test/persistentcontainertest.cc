// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief A unit test for the PersistentContainer
 */

#include <config.h>

#include <iostream>
#include <cassert>

#include <dune/common/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
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

  typedef typename GridType::LeafGridView GridView;
  const GridView &view = grid.leafView();
  typedef typename GridView::template Codim<0>::Iterator EIterator;

  {
    const EIterator &eend = view.template end<0>();
    for(EIterator eit = view.template begin<0>(); eit != eend; ++eit)
    {
      container0[*eit] = eit->geometry().center();
      const Dune::GenericReferenceElement< typename GridType::ctype, GridType::dimension > &refElement
        = Dune::GenericReferenceElements< typename GridType::ctype, GridType::dimension >::general( eit->type() );
      for (int i=0; i<eit->template count<1>(); ++i)
        container1(*eit,i) = eit->geometry().global( refElement.position(i,1) );
      for (int i=0; i<eit->template count<2>(); ++i)
        container2(*eit,i) = eit->geometry().global( refElement.position(i,2) );
    }
  }

  grid.globalRefine(1);
  container0.update();
  container1.update();
  container2.update();

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
      typename EIterator::EntityPointer::Entity::EntityPointer up = eit->father();
      while ( !container0[*up].used )
      {
        if (up->level() == 0)
        {
          std::cout << "ERROR: could not find a father element in container - stop testing" << std::endl;
          ret = false;
          break;
        }
        up = up->father();
      }
      if ( ( container0[*up].coord - up->geometry().center() ).two_norm() > 1e-8 )
      {
        std::cout << "ERROR: wrong data stored in container0 - stop testing" << std::endl;
        ret = false;
        break;
      }
      const Dune::GenericReferenceElement< typename GridType::ctype, GridType::dimension > &refElement
        = Dune::GenericReferenceElements< typename GridType::ctype, GridType::dimension >::general( eit->type() );
      for (int i=0; i<eit->template count<1>(); ++i)
        if ( ( container1(*up,i).coord - up->geometry().global( refElement.position(i,1) ) ).two_norm() > 1e-8 )
        {
          std::cout << "ERROR: wrong data stored in container1 - stop testing" << std::endl;
          ret = false;
          break;
        }
      for (int i=0; i<eit->template count<2>(); ++i)
        if ( ( container2(*up,i).coord - up->geometry().global( refElement.position(i,2) ) ).two_norm() > 1e-8 )
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
  MPIHelper & mpihelper = MPIHelper::instance(argc,argv);

  // /////////////////////////////////////////////////////////////////////////////
  //   Test YaspGrid
  // /////////////////////////////////////////////////////////////////////////////
  {
    typedef YaspGrid<2> GridType;
    Dune::FieldVector<double,2> Len; Len = 1.0;
    Dune::FieldVector<int,2> s; s = 2; s[0] = 6;
    Dune::FieldVector<bool,2> p; p = false;
    int overlap = 1;
    GridType grid(Len,s,p,overlap);
    std::cout << "Testing YaspGrid" << std::endl;
    test(grid);
  }

#if HAVE_ALUGRID
  {
    typedef ALUCubeGrid<2,2> GridType;
    array<unsigned int,2> elements2d;
    elements2d.fill(4);
    shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(FieldVector<double,2>(0),
                                                                                FieldVector<double,2>(1), elements2d);
    std::cout << "Testing ALUGrid" << std::endl;
    test(*grid);
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
