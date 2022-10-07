// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#define VERIFY(t,msg) do { if (!((t))) DUNE_THROW(Dune::Exception, "Check " #t " failed (" msg ")"); } while (false)

template<typename GV, typename R>
void checkRange(GV gv, int codim, R range, bool verify_less = false)
{
  int count = 0;
  int indices = 0;
  for (auto& e : range)
    {
      VERIFY(codim == e.codimension,"wrong codimension");
      ++count;
      indices += gv.indexSet().index(e);
    }
  if (verify_less)
    VERIFY(count <= gv.size(codim),"wrong number of elements for codim " << codim << " ");
  else
    VERIFY(count == gv.size(codim),"wrong number of elements for codim " << codim << " ");
}


// The grid *must* be refined at least once!
template<typename Grid, typename OS>
void check_ranges(const Grid& grid, OS&& os) {

  const int dim = Grid::dimension;

  auto gv = grid.leafGridView();

  os << "Checking entity ranges with default PartitionSet... ";

  // check elements
  checkRange(gv,0,elements(gv));
  // check facets
  checkRange(gv,1,facets(gv));
  // check edges
  checkRange(gv,dim-1,edges(gv));
  // check vertices
  checkRange(gv,dim,vertices(gv));
  // check elements using codim
  checkRange(gv,0,entities(gv,Dune::Codim<0>()));
  // check facets using codim
  checkRange(gv,1,entities(gv,Dune::Codim<1>()));
  // check edges using dim
  checkRange(gv,dim-1,entities(gv,Dune::Dim<1>()));
  // check vertices using dim
  checkRange(gv,dim,entities(gv,Dune::Dim<0>()));

  os << "OK" << std::endl;
  os << "Checking entity ranges with non-default PartitionSet... ";

  // check versions with PartitionSet parameter
  checkRange(gv,0,elements(gv,Dune::Partitions::interiorBorder),true);
  checkRange(gv,1,facets(gv,Dune::Partitions::all),true);
  checkRange(gv,dim-1,edges(gv,Dune::Partitions::interior + Dune::Partitions::border),true);
  checkRange(gv,dim,vertices(gv,Dune::Partitions::interiorBorder + Dune::Partitions::overlap + Dune::Partitions::front),false);
  checkRange(gv,1,entities(gv,Dune::Codim<1>(),Dune::Partitions::interiorBorder),true);
  checkRange(gv,dim-1,entities(gv,Dune::Dim<1>(),Dune::Partitions::all - Dune::Partitions::ghost),true);

  os << "OK" << std::endl;
  os << "Checking hierarchic entity range... ";

  // check hierarchic entity range
  {
    int count = 0;
    int indices = 0;

    // we need a LevelGridView higher up in the hierarchy for this
    auto gv = grid.levelGridView(grid.maxLevel()-1);

    for (auto& e : descendantElements(*elements(gv).begin(),grid.maxLevel()))
      {
        ++count;
        indices += grid.levelGridView(e.level()).indexSet().index(e);
      }

    VERIFY(count == 1 << dim,"wrong number of descendant elements");
  }

  os << "OK" << std::endl;
  os << "Checking intersection range... ";

  // check intersections
  {
    int count = 0;
    int indices = 0;

    for (auto& is : intersections(gv,*elements(gv).begin()))
      {
        ++count;
        indices += is.indexInInside();
      }

    VERIFY(count == 2 * dim,"wrong number of intersections");
  }

  os << "OK" << std::endl;
}


// little helper to only print output on rank 0

template<typename S, typename C>
struct Rank0Stream
{

  Rank0Stream(S& s, C c)
    : _s(s)
    , _c(c)
  {}

  template<typename T>
  Rank0Stream& operator<<(const T& t)
  {
    if (_c.rank() == 0)
      _s << t;
    return *this;
  }

  // needed for manipulators (std::endl etc.)
  Rank0Stream& operator<<(std::ostream& (*f)(std::ostream&))
  {
    if (_c.rank() == 0)
      f(_s);
    return *this;
  }

  S& _s;
  C _c;

};

template<typename S, typename C>
Rank0Stream<S,C> rank0Stream(S& s, C c)
{
  return {s,c};
}

template<typename OS>
void check_yasp_3d(OS&& os)
{

  os << "Running tests on 3D YaspGrid..." << std::endl;

  const int dim = 3;

  Dune::FieldVector<double,dim> Len(1.0);
  std::array<int,dim> s;
  std::fill(s.begin(), s.end(), 8);
  std::bitset<dim> p;
  int overlap = 1;

  Dune::YaspGrid<dim> grid(Len,s,p,overlap);

#if HAVE_MPI
  os << "Parallel run on " << grid.comm().size() << " processors" << std::endl;
#else
  os << "Sequential run" << std::endl;
#endif

  // refine once so that we can check the hierarchic iterator range
  grid.globalRefine(1);

  check_ranges(grid,os);
}


int main(int argc , char **argv) {
  try {
    // Initialize MPI, if present

    Dune::MPIHelper::instance(argc, argv);
    check_yasp_3d(rank0Stream(std::cout,Dune::MPIHelper::getCommunication()));

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
