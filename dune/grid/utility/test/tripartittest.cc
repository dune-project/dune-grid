// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the VertexOrder classes
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/albertagrid.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/tripartit.hh>
#include <dune/grid/yaspgrid.hh>

void fail(int &result) {
  result = 1;
}
void pass(int &result) {
  if(result == 77) result = 0;
}

template<class GV, class Partitioning>
class PartitioningVTKAdaptor :
  public Dune::VTKFunction<GV>
{
  typedef Dune::VTKFunction<GV> Base;

  const Partitioning &partitioning_;

public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  enum { dim = Base::dim };

  PartitioningVTKAdaptor(const Partitioning &partitioning) :
    partitioning_(partitioning)
  { }

  int ncomps() const
  {
    return 1;
  }

  double evaluate(int comp, const Entity &e,
                  const Dune::FieldVector< ctype, dim > &xi) const
  {
    return partitioning_.getPartition(e);
  }
  std::string name() const
  {
    return "partition";
  }
};

template<class GV, class Partitioner>
class ColoringVTKAdaptor :
  public Dune::VTKFunction<GV>
{
  typedef Dune::VTKFunction<GV> Base;
  typedef typename Base::ctype ctype;
  enum { dim = Base::dim };

  const Partitioner &partitioner_;

public:
  typedef typename Base::Entity Entity;

  ColoringVTKAdaptor(const Partitioner &partitioner) :
    partitioner_(partitioner)
  { }

  int ncomps() const
  {
    return 1;
  }

  double evaluate(int comp, const Entity &e,
                  const Dune::FieldVector< ctype, dim > &xi) const
  {
    return partitioner_.color(partitioner_.partitioning().getPartition(e));
  }
  std::string name() const
  {
    return "color";
  }
};

template<class GV>
void testTripartitColoring(const GV &gv, std::size_t overlap, int &result,
                           const std::string &vtkPrefix = "")
{
  pass(result);

  typedef Dune::GeneralFilteredPartitioning<GV> Partitioning;
  Partitioning partitioning(gv);

  typedef Dune::RecursiveEquidistantPartitioner<GV, Partitioning> Partitioner;
  Partitioner partitioner(gv, partitioning, overlap);

  while(partitioner.globalRefine())
    std::cout << "Number of Partitions: " << partitioning.partitions()
              << std::endl;

  if(vtkPrefix != "")
  {
    Dune::VTKWriter<GV> writer(gv);
    writer.addCellData
      (new PartitioningVTKAdaptor<GV, Partitioning>(partitioning));
    writer.addCellData(new ColoringVTKAdaptor<GV, Partitioner>(partitioner));
    writer.write(vtkPrefix);
  }
}

template<int dim>
void testYasp(int &result) {
  typedef Dune::YaspGrid<dim> Grid;

  Dune::array<int, dim> s;
  std::fill(s.begin(), s.end(), 1);

  Grid grid(Dune::FieldVector<typename Grid::ctype, dim>(1),
            s, std::bitset<dim>(), 0);
  grid.globalRefine(8);
  testTripartitColoring(grid.leafView(), 1, result, "tripartit-yasp");
}

std::string removeSuffix(const std::string &str, const std::string &suffix)
{
  if(str.size() > suffix.size() &&
     str.substr(str.size()-suffix.size()) == suffix)
    return str.substr(0, str.size()-suffix.size());
  else
    return str;
}

#if HAVE_ALUGRID || HAVE_ALBERTA
template<int dim>
void testGmsh(const std::string &fname, int &result) {
#if HAVE_ALUGRID
  typedef Dune::ALUSimplexGrid<dim, dim> Grid;
#else // HAVE_ALBERTA
  typedef Dune::AlbertaGrid<dim, dim> Grid;
#endif

  Dune::shared_ptr<Grid> gridp(Dune::GmshReader<Grid>::read(fname));
  testTripartitColoring(gridp->leafView(), 1, result,
                        removeSuffix(fname, ".msh"));
}
#endif // HAVE_ALUGRID || HAVE_ALBERTA

int main (int argc , char **argv)
try {

  // this method calls MPI_Init, if MPI is enabled
  Dune::MPIHelper::instance(argc,argv);

  int result = 77;

  if(argc > 1)
  {
#if HAVE_ALUGRID || HAVE_ALBERTA
    testGmsh<2>(argv[1], result);
#endif // HAVE_ALUGRID
  }
  else
    testYasp<2>(result);

  return result;
}
catch (const Dune::Exception &e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
