// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/gnuplot.hh>

template <class G, class IS>
void testIO(const G & g, const IS & is, std::string fname)
{
  // create data buffer
  std::vector<double> celldata;
  std::vector<float> vertexdata;

  for (int i=0; i<is.size(0); i++)
  {
    celldata.push_back(1+ 0.22*i);
    vertexdata.push_back(0.22*i);
  }
  vertexdata.push_back(0.22*is.size(0));

  // create writer
  Dune::GnuplotWriter<G,IS> gnuplot(g,is);

  // register data
  gnuplot.template addVertexData< std::vector<float> >(vertexdata, std::string("vertexdata"));
  gnuplot.template addCellData(celldata, std::string("celldata"));

  // write data
  gnuplot.write(fname);
}

int main()
{
  try {
    int n[] = { 10 };
    double h[] = { 1.0 };
    Dune::SGrid<1,1> sgrid(n,h);

    testIO(sgrid, sgrid.leafIndexSet(), "sgrid-leaf.data");
    testIO(sgrid, sgrid.levelIndexSet(0), "sgrid-level.data");
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  return 0;
}
