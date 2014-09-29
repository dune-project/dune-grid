// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

// make it possible to compile fallback code on newer compilers
#ifdef CHECK_WITHOUT_RANGE_BASED_FOR
#undef HAVE_RANGE_BASED_FOR
#endif

#define VERIFY(t,msg) do { if (!((t))) DUNE_THROW(Dune::Exception, "Check " #t " failed (" msg ")"); } while (false)

#if HAVE_RANGE_BASED_FOR

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

#else // HAVE_RANGE_BASED_FOR

template<typename GV, typename R>
void checkRange(GV gv, int codim, R range, bool verify_less = false)
{
  int count = 0;
  int indices = 0;
  for (typename R::iterator it = range.begin(); it != range.end(); ++it)
    {
      VERIFY(codim == it->codimension,"wrong codimension");
      ++count;
      indices += gv.indexSet().index(*it);
    }
  if (verify_less)
    VERIFY(count <= gv.size(codim),"wrong number of elements for codim " << codim << " ");
  else
    VERIFY(count == gv.size(codim),"wrong number of elements for codim " << codim << " ");
}

template<typename R>
void checkIntersectionRange(const R& range, int& count, int& indices)
{
  for (typename R::iterator it = range.begin(); it != range.end(); ++it)
    {
      ++count;
      indices += it->indexInInside();
    }
}

template<typename R, typename Grid>
void checkDescendantElementRange(const R& range, const Grid& grid, int& count, int& indices)
{
  for (typename R::iterator it = range.begin(); it != range.end(); ++it)
    {
      ++count;
      indices += grid.levelGridView(it->level()).indexSet().index(*it);
    }
}

#endif // HAVE_RANGE_BASED_FOR

template<typename OS>
void check_yasp(OS&& os) {

  const int dim = 3;

  Dune::FieldVector<double,dim> Len(1.0);
  Dune::array<int,dim> s;
  std::fill(s.begin(), s.end(), 8);
  std::bitset<dim> p;
  int overlap = 1;

#if HAVE_MPI
  Dune::YaspGrid<dim> grid(MPI_COMM_WORLD,Len,s,p,overlap);
  os << "Parallel run on " << grid.comm().size() << " processors" << std::endl;
#else
  Dune::YaspGrid<dim> grid(Len,s,p,overlap);
  os << "Sequential run" << std::endl;
#endif

  // refine once so that we can check the hierarchic iterator range
  grid.globalRefine(1);

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

#if HAVE_RANGE_BASED_FOR

    for (auto& e : descendantElements(*elements(gv).begin(),grid.maxLevel()))
      {
        ++count;
        indices += grid.levelGridView(e.level()).indexSet().index(e);
      }

#else // HAVE_RANGE_BASED_FOR

    // do the check in a function to capture the return type of the descendantElements call
    // -> greatly simplifies code
    checkDescendantElementRange(descendantElements(*(elements(gv).begin()),grid.maxLevel()),grid,count,indices);

#endif // HAVE_RANGE_BASED_FOR

    VERIFY(count == 1 << dim,"wrong number of descendant elements");
  }

  os << "OK" << std::endl;
  os << "Checking intersection range... ";

  // check intersections
  {
    int count = 0;
    int indices = 0;

#if HAVE_RANGE_BASED_FOR

    for (auto& is : intersections(gv,*elements(gv).begin()))
      {
        ++count;
        indices += is.indexInInside();
      }

#else // HAVE_RANGE_BASED_FOR

    // do the check in a function to capture the return type of the intersections call
    // -> greatly simplifies code
    checkIntersectionRange(intersections(gv,*(elements(gv).begin())),count,indices);

#endif // HAVE_RANGE_BASED_FOR

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



int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);
    check_yasp(rank0Stream(std::cout,Dune::MPIHelper::getCollectiveCommunication()));

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
