// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_TEST_YASPGRID_HH
#define DUNE_GRID_TEST_TEST_YASPGRID_HH

#include <string>
#include <memory>

#include <dune/grid/yaspgrid.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkcommunicate.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkadaptation.hh>
#include <dune/grid/test/checkpartition.hh>

#include <dune/grid/test/checkcomcorrectness.hh>

#include <dune/grid/yaspgrid/coordinates.hh>


template<int dim, class CC>
struct YaspFactory
{};

template<int dim>
struct YaspFactory<dim, Dune::EquidistantCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim>* buildGrid(
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false,
      bool useGenericConstructor = false)
  {
    std::cout << " using equidistant coordinate container";
    if (useGenericConstructor)
      std::cout << " with generic constructor";
    std::cout << "!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> len(1.0);
    std::array<int,dim> s;
    std::fill(s.begin(), s.end(), 4);
    s[0] = 8;
    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    Dune::YaspGrid<dim>* grid;

    if (useGenericConstructor)
    {
      Dune::EquidistantCoordinates<double,dim> coordinates(len,s);
      grid = new Dune::YaspGrid<dim>(coordinates,p,overlap);
    }
    else
      grid = new Dune::YaspGrid<dim>(len,s,p,overlap);

    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template<int dim>
struct YaspFactory<dim, Dune::EquidistantOffsetCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >* buildGrid(
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false,
      bool useGenericConstructor = false)
  {
    if (useGenericConstructor)
      std::cout << " using equidistant coordinate container with non-zero origin and generic constructor!" << std::endl << std::endl;
    else
      std::cout << " using equidistant coordinate container with non-zero origin!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> lowerleft(-1.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    std::array<int,dim> s;
    std::fill(s.begin(), s.end(), 4);
    s[0] = 8;
    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >* grid;

    if (useGenericConstructor)
    {
      Dune::EquidistantOffsetCoordinates<double,dim> coordinates(lowerleft,upperright,s);
      grid = new Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >(coordinates,p,overlap);
    }
    else
      grid = new Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >(lowerleft,upperright,s,p,overlap);

    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template<int dim>
struct YaspFactory<dim, Dune::TensorProductCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >* buildGrid(
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false,
      bool useGenericConstructor = false)
  {
    if (useGenericConstructor)
      std::cout << " using tensorproduct coordinate container and generic constructor!" << std::endl << std::endl;
    else
      std::cout << " using tensorproduct coordinate container!" << std::endl << std::endl;

    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    std::array<std::vector<double>,dim> coords;
    coords[0].resize(9);
    coords[0][0] = -1.0;
    coords[0][1] = -0.5;
    coords[0][2] = -0.25;
    coords[0][3] = -0.125;
    coords[0][4] =  0.0;
    coords[0][5] =  0.125;
    coords[0][6] =  0.25;
    coords[0][7] =  0.5;
    coords[0][8] =  1.0;
    for (int i=1; i<dim; i++)
    {
      coords[i].resize(5);
      coords[i][0] = -1.0;
      coords[i][1] = -0.5;
      coords[i][2] = 0.0;
      coords[i][3] = 0.5;
      coords[i][4] =  1.0;
    }

    Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >* grid;

    if (useGenericConstructor)
    {
      std::array<int,dim> offset;
      std::fill(offset.begin(), offset.end(), 0);
      Dune::TensorProductCoordinates<double,dim> coordinates(coords,offset);
      grid = new Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >(coordinates,p,overlap);
    }
    else
      grid = new Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >(coords,p,overlap);

    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template <int dim, class CC>
void check_yasp(std::string testID, Dune::YaspGrid<dim,CC>* grid) {
  std::cout << std::endl << "YaspGrid<" << dim << ">";

  gridcheck(*grid);

  checkIterators ( grid->leafGridView() );
  checkIterators ( grid->levelGridView(0) );

  // check communication interface
  checkCommunication(*grid,-1,Dune::dvverb);
  for(int l=0; l<=grid->maxLevel(); ++l)
    checkCommunication(*grid,l,Dune::dvverb);

  // check communication correctness
  Dune::GridCheck::check_communication_correctness(grid->leafGridView());

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
  filename << testID << "-output" <<grid->comm().rank();
  file.open(filename.str());
  file << *grid << std::endl;
  file.close();

  delete grid;
}

template <int dim, class CC = Dune::EquidistantCoordinates<double,dim> >
void check_backuprestore(std::string testID, Dune::YaspGrid<dim,CC>* grid)
{
   typedef Dune::YaspGrid<dim,CC> Grid;
   grid->globalRefine(2);

   Dune::BackupRestoreFacility<Grid>::backup(*grid, testID+"-backup");

   // avoid that processes that having nothing to backup try to restore
   // a grid that has not been backuped yet.
   grid->comm().barrier();
   Grid* restored = Dune::BackupRestoreFacility<Grid>::restore(testID+"-backup");

   // write a backup of the restored file. this has to be identical to backup
   Dune::BackupRestoreFacility<Grid>::backup(*restored, testID+"-copy");

   if ((std::is_same<CC,Dune::TensorProductCoordinates<double,dim> >::value) || (grid->comm().rank() == 0))
   {
     // check whether copy and backup are equal
     std::ostringstream s1,s2;
     s1 << testID << "-backup";
     s2 << testID << "-copy";
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

   check_yasp(testID, restored);

   delete grid;
}

#endif
