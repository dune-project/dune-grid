// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERUG_HH
#define DUNE_DGFPARSERUG_HH
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl<UGGrid<dim,dimworld> > {
  public:
    static UGGrid<dim,dimworld>* generate(MacroGrid& mg,
                                          const char* filename,int MPICOMM=-1);
  };
}
#include "dgfug.cc"
#endif
#endif
