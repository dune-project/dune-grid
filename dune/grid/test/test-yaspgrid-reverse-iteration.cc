// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

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
    std::array<int,dim> s;
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
    std::array<int,dim> s;
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

    std::array<std::vector<double>,dim> coords;
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
  std::cout << std::endl << "YaspGrid<" << dim << ">" << std::endl;

  const int codim = 1;

  if (grid == NULL)
    grid = YaspFactory<dim,CC>::buildGrid();

  auto gv = grid->leafGridView();

  auto it = gv.template begin<codim>();
  auto end = gv.template end<codim>();
  auto rit = gv.template rbegin<codim>();
  auto rend = gv.template rend<codim>();
  auto& is = gv.indexSet();

  std::vector<typename decltype(gv)::template Codim<codim>::Entity> entities_;

  for ( ; it != end ; ++it)
    entities_.push_back(*it);

  auto hit = entities_.rbegin();

  for(; rit != rend ; ++rit, ++hit)
    {
      std::cout << std::setw(4) << is.index(*hit) << " " << std::setw(4) << is.index(*rit) << std::endl;
      if (hit->geometry().center() != rit->geometry().center())
        std::cout << "Error: " << hit->geometry().center() << " != " << rit->geometry().center() << std::endl;
    }

  for(const auto& e : facets(gv,Dune::Direction::reverse))
    std::cout << is.index(e) << " ";
  std::cout << std::endl;

  delete grid;
}


int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);

    check_yasp(YaspFactory<3,Dune::EquidistantOffsetCoordinates<double,3> >::buildGrid());

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
