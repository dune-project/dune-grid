// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERS_HH
#define DUNE_DGFPARSERS_HH
#include <dune/grid/sgrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl< SGrid<dim,dimworld> > {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static SGrid<dim,dimworld>* generate(MacroGrid& mg,
                                         const char* filename,
                                         MPICommunicatorType MPICOMM = MPIHelper::getCommunicator());
  };
  template <int dim, int dimworld>
  struct DGFGridInfo< SGrid<dim,dimworld> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dim);}
  };
}
#include "dgfs.cc"
#endif
