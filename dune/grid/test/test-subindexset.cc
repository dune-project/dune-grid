#include <config.h>

#include <dune/grid/common/subindexset.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
// #include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk.hh>

#include <iostream>

using namespace Dune;

template<class GV>
bool test_grid_view(const GV& gv)
{
  return test_index_set(gv,gv.indexSet());
}


int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    const auto &helper = Dune::MPIHelper::instance(argc, argv);

    std::array<unsigned int,2> elements2d;
    elements2d.fill(3);

    FieldVector<double,2> lowerLeft(0);
    FieldVector<double,2> upperRight(1);

    using Grid = Dune::YaspGrid<2>;

    // Test creation of 2d uggrid
    auto grid_ptr = StructuredGridFactory<Grid>::createCubeGrid(lowerLeft,upperRight,elements2d);

    Dune::gridinfo(*grid_ptr);
    auto gv = grid_ptr->leafGridView();
    using GV = decltype(gv);

    Dune::VTKWriter<GV> writer(gv);
    writer.write("anything",VTK::OutputType::ascii);

    Dune::PDELab::EntitySetIndexSet<GV> is(gv);

    for (auto i : is._gt_off_base)
      std::cout << i << " ";
    std::cout << std::endl;

    // print all the index set
    for(auto&& vertices : vertices(gv))
      std::cout << is.unique_base_index(vertices) << " ";
    for(auto&& edges : edges(gv))
      std::cout << is.unique_base_index(edges) << " ";
    for(auto&& element : elements(gv))
      std::cout << is.unique_base_index(element) << " ";
    std::cout << std::endl;

    auto vertex = gv.template begin<2>();
    // is.unmark(0);
    // is.unmark(1);
    // is.unmark(2);
    // is.unmark(*vertex);
    is.unmark(*std::next(vertex));
    is.unmark(*std::next(std::next(vertex)));
    std::cout << is._marked_gt.to_string() << std::endl;
    std::cout << is._active_gt.to_string() << std::endl;
    // assert(is.getMark(*vertex) == true);
    // assert(is.getMark(*std::next(vertex)) == false);
    // assert(is.getMark(*std::next(std::next(vertex))) == true);
    auto edge = gv.template begin<1>();
    // assert(is.getMark(*edge) == false);
    auto element = gv.template begin<0>();
    // assert(is.getMark(*element) == false);

    is.update(true);

    for (auto i : is._gt_off)
      std::cout << i << " ";
    std::cout << std::endl;

    for (auto [key_i,value_i] : is._index_off)
        for (auto [key_j,value_j] : value_i)
      std::cout << "["<< key_j << "," << value_j << "] ";
    std::cout << std::endl;

    std::cout << is.contains(*vertex) << std::endl;
    std::cout << is.contains(*std::next(vertex)) << std::endl;
    std::cout << is.contains(*std::next(std::next(vertex))) << std::endl;
    std::cout << is.contains(*std::next(std::next(std::next(vertex)))) << std::endl;

    std::cout << is.unique_index(*vertex) << std::endl;
    std::cout << is.unique_index(*std::next(vertex)) << std::endl;
    std::cout << is.unique_index(*std::next(std::next(vertex))) << std::endl;
    std::cout << is.unique_index(*std::next(std::next(std::next(vertex)))) << std::endl;



    // is.update(true,Dune::Partitions::interior);


    // is.update(true,Dune::Partitions::overlap);

    // for (auto [key,value] : is._index_off)
    //   std::cout << "["<< key << "," << value << "] ";
    // std::cout << std::endl;
    //
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}