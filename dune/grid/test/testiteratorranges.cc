#include "config.h"

#include <dune/grid/sgrid.hh>

#define VERIFY(t,msg) do { if (!((t))) DUNE_THROW(Dune::Exception, "Check " #t " failed (" msg ")"); } while (false)

#if HAVE_RANGE_BASED_FOR

template<typename GV, typename R>
void checkRange(GV gv, int codim, R range)
{
  int count = 0;
  int indices = 0;
  for (auto& e : range)
    {
      ++count;
      indices += gv.indexSet().index(e);
    }
  VERIFY(count == gv.size(codim),"wrong number of elements at codim " << codim << " ");
}

#else // HAVE_RANGE_BASED_FOR

template<typename GV, typename R>
void checkRange(GV gv, int codim, R range)
{
  int count = 0;
  int indices = 0;
  for (typename R::iterator it = range.begin(); it != range.end(); ++it)
    {
      ++count;
      indices += gv.indexSet().index(*it);
    }
  VERIFY(count == gv.size(codim),"wrong number of elements at codim " << codim << " ");
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

template<typename R>
void checkChildElementRange(const R& range, int& count, int& indices)
{
  for (typename R::iterator it = range.begin(); it != range.end(); ++it)
    {
      ++count;
      indices += it->level();
    }
}

#endif // HAVE_RANGE_BASED_FOR

int main(int argc, char** argv)
{
  try {

    int n[] = { 1, 1 ,1 };
    double h[] = { 1.0, 1.0, 1.0 };

    typedef Dune::SGrid<3,3> Grid;
    Grid grid(n,h);

    grid.globalRefine(4);

    Grid::LeafGridView gv(grid.leafGridView());

    elements(gv);

    // check elements
    checkRange(gv,0,elements(gv));
    // check faces
    checkRange(gv,1,faces(gv));
    // check edges
    checkRange(gv,2,edges(gv));
    // check vertices
    checkRange(gv,3,vertices(gv));
    // check elements using codim
    checkRange(gv,0,entities(gv,Dune::Codim<0>()));
    // check elements using dim
    checkRange(gv,1,entities(gv,Dune::Dim<2>()));
    // check vertices using codim
    checkRange(gv,2,entities(gv,Dune::Codim<2>()));
    // check vertices using dim
    checkRange(gv,2,entities(gv,Dune::Dim<1>()));

    // check intersections
    {
      int count = 0;
      int indices = 0;

#if HAVE_RANGE_BASED_FOR

      for (auto& is : intersections(*elements(gv).begin(),gv))
        {
          ++count;
          indices += is.indexInInside();
        }

#else // HAVE_RANGE_BASED_FOR

      // do the check in a function to capture the return type of the intersections call
      // -> greatly simplifies code
      checkIntersectionRange(intersections(*(elements(gv).begin()),gv),count,indices);

#endif // HAVE_RANGE_BASED_FOR

      VERIFY(count == 4,"wrong number of intersections");
    }

    // check hierarchic iterators
    {
      int count = 0;
      int indices = 0;

#if HAVE_RANGE_BASED_FOR

      for (auto& e : childElements(*elements(grid.levelGridView(0)).begin(),3))
        {
          ++count;
          indices += e.level();
        }

#else // HAVE_RANGE_BASED_FOR

      // do the check in a function to capture the return type of the childElements call
      // -> greatly simplifies code
      checkChildElementRange(childElements(*elements(grid.levelGridView(0)).begin(),3),count,indices);

#endif // HAVE_RANGE_BASED_FOR

      VERIFY(count > 0,"wrong number of child elements");
    }

  } catch (Dune::Exception& e) {
    std::cout << "Dune exception: " << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }
  return 0;
}
