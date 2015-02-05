// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/backuprestore.hh>
#include <dune/grid/yaspgrid/factory.hh>

#include "gridcheck.cc"
#include "checkcommunicate.hh"
#include "checkgeometryinfather.hh"
#include "checkintersectionit.hh"
#include "checkiterators.hh"
#include "checkadaptation.hh"
#include "checkpartition.cc"

template<int dim, class CC>
struct YaspFactory
{};

template<int dim>
struct YaspFactory<dim, Dune::EquidistantCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim>* buildGrid()
  {
    std::cout << " using equidistant coordinate container!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> Len(1.0);
    Dune::array<int,dim> s;
    std::fill(s.begin(), s.end(), 8);
    std::bitset<dim> p(0);
    int overlap = 1;

    return new Dune::YaspGrid<dim>(Len,s,p,overlap);
  }
};

template<int dim>
struct YaspFactory<dim, Dune::EquidistantOffsetCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >* buildGrid()
  {
    std::cout << " using equidistant coordinate container with non-zero origin!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> lowerleft(-1.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    Dune::array<int,dim> s;
    std::fill(s.begin(), s.end(), 8);
    std::bitset<dim> p(0);
    int overlap = 1;

    return new Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >(lowerleft,upperright,s,p,overlap);
  }
};

template<int dim>
struct YaspFactory<dim, Dune::TensorProductCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >* buildGrid()
  {
    std::cout << " using tensorproduct coordinate container!" << std::endl << std::endl;

    std::bitset<dim> p(0);
    int overlap = 1;

    Dune::array<std::vector<double>,dim> coords;
    for (int i=0; i<dim; i++)
    {
      coords[i].resize(9);
      coords[i][0] = -1.0;
      coords[i][1] = -0.5;
      coords[i][2] = -0.25;
      coords[i][3] = -0.125;
      coords[i][4] =  0.0;
      coords[i][5] =  0.125;
      coords[i][6] =  0.25;
      coords[i][7] =  0.5;
      coords[i][8] =  1.0;
    }

    return new Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >(coords,p,overlap);
  }
};

template <int dim, class CC>
void check_yasp(Dune::YaspGrid<dim,CC>* grid) {
  std::cout << std::endl << "YaspGrid<" << dim << ">";

  if (grid == NULL)
    grid = YaspFactory<dim,CC>::buildGrid();

  gridcheck(*grid);
  //grid->globalRefine(2);

  gridcheck(*grid);

  checkIterators ( grid->leafGridView() );
  checkIterators ( grid->levelGridView(0) );

  // check communication interface
  checkCommunication(*grid,-1,Dune::dvverb);
  for(int l=0; l<=grid->maxLevel(); ++l)
    checkCommunication(*grid,l,Dune::dvverb);

  // check geometry lifetime
  checkGeometryLifetime( grid->leafGridView() );
  // check the method geometryInFather()
  checkGeometryInFather(*grid);
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(*grid);
  // check grid adaptation interface
  checkAdaptRefinement(*grid);
  checkPartitionType( grid->leafGridView() );

  std::ofstream file;
  std::ostringstream filename;
  filename << "output" <<grid->comm().rank();
  file.open(filename.str());
  file << *grid << std::endl;
  file.close();

  delete grid;
}

template <int dim, class CC = Dune::EquidistantCoordinates<double,dim> >
void check_backuprestore(Dune::YaspGrid<dim,CC>* grid)
{
   typedef Dune::YaspGrid<dim,CC> Grid;
   grid->globalRefine(2);

   Dune::BackupRestoreFacility<Grid>::backup(*grid, "backup");

   // avoid that processes that having nothing to backup try to restore
   // a grid that has not been backuped yet.
   grid->comm().barrier();
   Grid* restored = Dune::BackupRestoreFacility<Grid>::restore("backup");

   // write a backup of the restored file. this has to be identical to backup
   Dune::BackupRestoreFacility<Grid>::backup(*restored, "copy");

   if ((std::is_same<CC,Dune::TensorProductCoordinates<double,dim> >::value) || (grid->comm().rank() == 0))
   {
     // check whether copy and backup are equal
     std::ostringstream s1,s2;
     s1 << "backup";
     s2 << "copy";
     if (std::is_same<CC,Dune::TensorProductCoordinates<double,dim> >::value)
     {
       s1 << grid->comm().rank();
       s2 << grid->comm().rank();
     }
     std::ifstream file1, file2;
     file1.open(s1.str());
     file2.open(s2.str());

     std::string token1, token2;
     while(!file1.eof() && !file2.eof())
     {
       file1 >> token1;
       file2 >> token2;
       if (token1 != token2)
         DUNE_THROW(Dune::Exception, "Error in BackupRestoreFacility");
     }
   }

   check_yasp(restored);

   delete grid;
}

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);

    check_yasp(YaspFactory<1,Dune::EquidistantCoordinates<double,1> >::buildGrid());
    check_yasp(YaspFactory<1,Dune::EquidistantOffsetCoordinates<double,1> >::buildGrid());
    check_yasp(YaspFactory<1,Dune::TensorProductCoordinates<double,1> >::buildGrid());

    check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid());
    check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid());
    check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid());

    check_yasp(YaspFactory<3,Dune::EquidistantCoordinates<double,3> >::buildGrid());
    check_yasp(YaspFactory<3,Dune::EquidistantOffsetCoordinates<double,3> >::buildGrid());
    check_yasp(YaspFactory<3,Dune::TensorProductCoordinates<double,3> >::buildGrid());

    // check the factory class for tensorproduct grids
    Dune::TensorYaspGridFactory<double,2> factory;
    factory.setStart(0,-100.);
    factory.fillIntervals(0,10,20.);
    factory.fillRange(0, 5, 130.);
    factory.geometricFillIntervals(0, 5, 2.0);

    factory.geometricFillRange(1,10,100.,1.,false);
    factory.fillRange(1,10,200);
    factory.geometricFillRange(1,10,250.,1.,true);
    factory.fillUntil(1,50,1000.);

    auto grid = factory.createGrid();
    delete grid;

    // check the backup restore facility
    check_backuprestore(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid());
    check_backuprestore(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid());
    check_backuprestore(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid());

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
