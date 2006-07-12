// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSERS_HH
#define DUNE_MACROGRIDPARSERS_HH
#include <dune/grid/sgrid.hh>
#include "macrogridparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl<SGrid<dim,dimworld> > {
  public:
    static SGrid<dim,dimworld>* generate(MacroGrid& mg,
                                         const char* filename,MPI_Comm MPICOMM=MPI_COMM_WORLD);
  };
}
#include "macrogrids.cc"
#endif
