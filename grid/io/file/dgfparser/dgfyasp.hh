// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH
#include <dune/grid/yaspgrid.hh>
#include "dgfparserblocks.hh"
#include "dgfparser.hh"
namespace Dune {
  template <int dim>
  class MacroGrid::Impl<YaspGrid<dim> > {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static YaspGrid<dim>* generate(MacroGrid& mg,
                                   const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };
  template <int dim>
  struct DGFGridInfo< YaspGrid<dim> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dim);}
  };
}
#include "dgfyasp.cc"
#endif
