// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERONED_HH
#define DUNE_DGFPARSERONED_HH
#include <dune/grid/onedgrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl<OneDGrid<dim,dimworld> > {
  public:
    static OneDGrid<dim,dimworld>* generate(MacroGrid& mg,
                                            const char* filename,int MPICOMM=-1);
  };
}
#include "dgfoned.cc"
#endif
