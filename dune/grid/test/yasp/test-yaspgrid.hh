// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_TEST_YASPGRID_HH
#define DUNE_GRID_TEST_TEST_YASPGRID_HH

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
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false)
  {
    std::cout << " using equidistant coordinate container!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> Len(1.0);
    std::array<int,dim> s;
    if (dim < 3)
      std::fill(s.begin(), s.end(), 8);
    else
      std::fill(s.begin(), s.end(), 4);
    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    auto grid = new Dune::YaspGrid<dim>(Len,s,p,overlap);
    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template<int dim>
struct YaspFactory<dim, Dune::EquidistantOffsetCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >* buildGrid(
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false)
  {
    std::cout << " using equidistant coordinate container with non-zero origin!" << std::endl << std::endl;

    Dune::FieldVector<double,dim> lowerleft(-1.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    std::array<int,dim> s;
    if (dim < 3)
      std::fill(s.begin(), s.end(), 8);
    else
      std::fill(s.begin(), s.end(), 4);
    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    auto grid = new Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim> >(lowerleft,upperright,s,p,overlap);
    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template<int dim>
struct YaspFactory<dim, Dune::TensorProductCoordinates<double,dim> >
{
  static Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >* buildGrid(
      bool keepPhysicalOverlap = true, int refCount = 0, bool periodic = false)
  {
    std::cout << " using tensorproduct coordinate container!" << std::endl << std::endl;

    std::bitset<dim> p(0);
    p[0] = periodic;
    int overlap = 1;

    std::array<std::vector<double>,dim> coords;
    if (dim < 3) {
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
    } else {
      for (int i=0; i<dim; i++)
      {
        coords[i].resize(7);
        coords[i][0] = -1.0;
        coords[i][1] = -0.5;
        coords[i][2] = -0.25;
        coords[i][3] =  0.0;
        coords[i][4] =  0.25;
        coords[i][5] =  0.5;
        coords[i][6] =  1.0;
      }
    }

    auto grid = new Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double,dim> >(coords,p,overlap);
    grid->refineOptions (keepPhysicalOverlap);
    grid->globalRefine (refCount);
    return grid;
  }
};

template <int dim, class CC>
void check_yasp(Dune::YaspGrid<dim,CC>* grid) {
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

#endif
