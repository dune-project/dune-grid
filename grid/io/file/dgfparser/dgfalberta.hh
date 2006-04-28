// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERALBERTA_HH
#define DUNE_DGFPARSERALBERTA_HH
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include "dgfparser.hh"
namespace Dune {
  template <int dim,int dimworld>
  class MacroGrid::Impl<AlbertaGrid<dim,dimworld> > {
  public:
    static AlbertaGrid<dim,dimworld>* generate(MacroGrid& mg,
                                               const char* filename,int MPICOMM=-1);
  };
}
#include "dgfalberta.cc"
#endif
#endif
