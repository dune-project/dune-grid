// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id:$

#include "config.h" // autoconf defines, needed by the dune headers

// dune headers
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <vector>
#include <unistd.h>

template<class G, class IS>
void doWrite(G & g, IS & is, Dune::VTKOptions::DataMode dm)
{
  static int run=0;
  enum { dim = G::dimension };

  Dune::VTKWriter<G, IS> vtk(g,is,dm);
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);
  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  run++;
  char name[256];
  snprintf(name,256,"vtktest-%i-ascii", run);
  vtk.write(name);
  unlink(name);

  snprintf(name,256,"vtktest-%i-binary", run);
  vtk.write(name, Dune::VTKOptions::binaryappended);
  unlink(name);
}

template<int dim>
void vtkCheck(int* n, double* h)
{
  std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;
  Dune::SGrid<dim,dim> g(n, h);
  g.globalRefine(1);

  doWrite(g,g.leafIndexSet(),Dune::VTKOptions::conforming);
  doWrite(g,g.leafIndexSet(),Dune::VTKOptions::nonconforming);
  doWrite(g,g.levelIndexSet(0),Dune::VTKOptions::conforming);
  doWrite(g,g.levelIndexSet(0),Dune::VTKOptions::nonconforming);
  doWrite(g,g.levelIndexSet(g.maxLevel()),Dune::VTKOptions::conforming);
  doWrite(g,g.levelIndexSet(g.maxLevel()),Dune::VTKOptions::nonconforming);
}

int main(int argc, char **argv)
{
  try {

    int n[] = { 5, 5, 5, 5 };
    double h[] = { 1.0, 2.0, 3.0, 4.0 };

    vtkCheck<1>(n,h);
    vtkCheck<2>(n,h);
    vtkCheck<3>(n,h);

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
