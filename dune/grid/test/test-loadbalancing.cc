// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#if HAVE_PARMETIS
#include <parmetis.h>
#  ifdef PARMETIS_MAJOR_VERSION
#    include <dune/grid/utility/parmetisgridpartitioner.hh>
#  endif
#endif

using namespace Dune;

template<int dim>
struct Ball {
  double radius;
  Dune::FieldVector<double, dim> center;

  Ball(const Dune::FieldVector<double, dim>& c, const double& r) : radius(r), center(c) {}

  double distanceTo(const Dune::FieldVector<double, dim>& other) const {
    return std::abs((center - other).two_norm() - radius);
  }
};

// Write a set of default parameters into a ParameterTree.
// The user can override all of these from the command line.
ParameterTree defaultParameters(int dim)
{
  ParameterTree result;
  if (dim==2)
  {
    result["n"] = "8 25";
    result["lower"] = "0 0";
    result["upper"] = "0.012 0.03";

    //result["center"] = "0.012 0.01"
    result["center"] = "0.006 0.005";
    //result["r"] = "0.004";
    result["r"] = "0.002";

    result["stepDisplacement"] = "0 0.001";
  }
  else if (dim==3)
  {
    result["n"] = "8 2 2";
    result["lower"] = "0 0 0";
    result["upper"] = "0.012 0.03 0.012";

    //result["center"] = "0.012 0.01 0.012"
    result["center"] = "0.006 0.005 0.006";
    result["r"] = "0.002";

    result["stepDisplacement"] = "0 0.001 0";
  }

  result["epsilon"] = "0.0008";
  result["levels"] = "2";
  result["steps"] = "30";

  return result;
}

// Define some constants and types
const int dim = 2;

typedef UGGrid<dim> GridType;
typedef GridType::LeafGridView GV;


int main(int argc, char** argv)
{
#if ! (HAVE_PARMETIS && defined(PARMETIS_MAJOR_VERSION))
  // Skip test -- without ParMetis it doesn't do anything useful
  std::cerr << "This test requires ParMetis and will be skipped.\n"
            << "Note that the emulation layer provided by scotch is not sufficient.\n";
  return 77;
#else
  // Create MPIHelper instance
  MPIHelper& mpihelper = MPIHelper::instance(argc, argv);

  if (0 == mpihelper.rank())
    std::cout << "Using " << mpihelper.size() << " Processes." << std::endl;

  ParameterTree parameterSet = defaultParameters(dim);

  // Create ug grid from structured grid
  const std::array<unsigned, dim> n = parameterSet.get<std::array<unsigned, dim> >("n");

  const FieldVector<double, dim>
    lower = parameterSet.get<FieldVector<double, dim> >("lower"),
    upper = parameterSet.get<FieldVector<double, dim> >("upper");

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, n);

  const GV gv = grid->leafGridView();

  // Create ball
  const FieldVector<double, dim> center = parameterSet.get<FieldVector<double, dim> >("center");
  const double r = parameterSet.get<double>("r");

  Ball<dim> ball(center, r);

  // Refine
  const size_t steps = parameterSet.get<size_t>("steps");
  const FieldVector<double, dim> stepDisplacement = parameterSet.get<FieldVector<double, dim> >("stepDisplacement");

  const double epsilon = parameterSet.get<double>("epsilon");
  const int levels = parameterSet.get<int>("levels");


  // Create initial partitioning using ParMETIS
  std::vector<unsigned> part(ParMetisGridPartitioner<GV>::partition(gv, mpihelper));

  // Transfer partitioning from ParMETIS to our grid
  grid->loadBalance(part, 0);

  for (size_t s = 0; s < steps; ++s) {
    std::cout << "Step " << s << " on " << mpihelper.rank() << " ..." << std::endl;

    for (int k = 0; k < levels; ++k) {
      std::cout << "   Refining level " << k << " on " << mpihelper.rank() << " ..." << std::endl;

      // select elements that are close to the sphere for grid refinement
      for (const auto& element : elements(gv, Partitions::interior))
      {
        if (ball.distanceTo(element.geometry().center()) < epsilon)
          grid->mark(1, element);
      }

      // adapt grid
      grid->adapt();

      // clean up markers
      grid->postAdapt();
    }

    mpihelper.getCommunication().barrier();

    // Repartition
    float itr = 1000; // ratio of inter-processor communication time compared to data redistribution time
                       // high ~> minimize edge-cut and have smaller communication time during calculations
                       // low  ~> do not move elements around between processes too much and thous reduce communication time during redistribution

    part = ParMetisGridPartitioner<GV>::repartition(gv, mpihelper, itr);

    // Transfer partitioning from ParMETIS to our grid
    grid->loadBalance(part, 0);

    // Output grid
    const std::string baseOutName = "RefinedGrid_";

    VTKWriter<GV> vtkWriter(gv);
    std::vector<int> rankField(gv.size(0));
    std::fill(rankField.begin(), rankField.end(), grid->comm().rank());
    vtkWriter.addCellData(rankField,"rank");
    vtkWriter.write(baseOutName+std::to_string(s));

    // If this is not the last step, move sphere and coarsen grid
    if (s+1 < steps) {
      // Move sphere a little
      ball.center += stepDisplacement;

      // Coarsen everything
      for (int k = 0; k < levels; ++k) {
        std::cout << "   Coarsening level " << k << " on " << mpihelper.rank() << " ..." << std::endl;

        for (const auto& element : elements(gv, Partitions::interior))
          grid->mark(-1, element);

        // adapt grid
        grid->adapt();

        // clean up markers
        grid->postAdapt();
      }
    }
  }

  return 0;
#endif   // (HAVE_PARMETIS && defined(PARMETIS_MAJOR_VERSION))
}
