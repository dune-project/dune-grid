// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include <iostream>
#include <ostream>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gnuplot.hh>

template <class GV>
void testIO(const GV & gridView, std::string fname)
{
  // create data buffer
  std::vector<double> celldata;
  std::vector<float> vertexdata;

  for (size_t i=0; i<gridView.indexSet().size(0); i++)
  {
    celldata.push_back(1+ 0.22*i);
    vertexdata.push_back(0.22*i);
  }
  vertexdata.push_back(0.22*gridView.indexSet().size(0));

  // create writer
  Dune::GnuplotWriter<GV> gnuplot(gridView);

  // register data
  gnuplot.template addVertexData< std::vector<float> >(vertexdata, std::string("vertexdata"));
  gnuplot.template addCellData(celldata, std::string("celldata"));

  // write data
  gnuplot.write(fname);
}

int main(int argc, char** argv)
{
  try
  {
    Dune::MPIHelper::instance(argc, argv);
    const unsigned int dim = 1;
    Dune::FieldVector<double,dim> length(1.0);
    std::array<int,dim> elements;
    elements[0] = 10;
    Dune::YaspGrid<1> grid(length, elements);

    testIO(grid.leafGridView(), "grid-leaf.data");
    testIO(grid.levelGridView(0), "grid-level.data");
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  return 0;
}
