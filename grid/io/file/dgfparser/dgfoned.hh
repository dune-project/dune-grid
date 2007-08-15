// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERONED_HH
#define DUNE_DGFPARSERONED_HH
#include <dune/grid/onedgrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <>
  class MacroGrid::Impl<OneDGrid> {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static OneDGrid* generate(MacroGrid& mg,
                              const char* filename,
                              MPICommunicatorType MPICOMM = MPIHelper::getCommunicator());
  };
  template <>
  struct DGFGridInfo<OneDGrid> {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return 0.5;}
  };
}
#include "dgfoned.cc"
#endif
