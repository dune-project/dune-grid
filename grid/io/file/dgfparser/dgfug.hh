// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERUG_HH
#define DUNE_DGFPARSERUG_HH

// only include if UG is used
#if defined ENABLE_UG
#include <dune/grid/uggrid.hh>
#include "dgfparser.hh"
namespace Dune
{

  /** \cond */
  template <int dim>
  class MacroGrid::Impl<UGGrid<dim> > {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static UGGrid<dim>* generate(MacroGrid& mg,
                                 const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };
  template <int dimw>
  struct DGFGridInfo< UGGrid<dimw> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return -1.;}
  };
  /** \endcond */

}
#include "dgfug.cc"
#endif

#endif
