// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH
#include <dune/grid/yaspgrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl<YaspGrid<dim,dimworld> > {
  public:
    static YaspGrid<dim,dimworld>* generate(MacroGrid& mg,
                                            const char* filename,int MPICOMM=-1);
  };
}
#include "dgfyasp.cc"
#endif
