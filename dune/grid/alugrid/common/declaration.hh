// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_FORWARDDECLARATION
#define DUNE_ALUGRID_FORWARDDECLARATION

#include <dune/common/collectivecommunication.hh>
#if HAVE_MPI
#include <dune/common/mpicollectivecommunication.hh>
#endif

namespace Dune {

  //! \brief basic element types for ALUGrid
  enum ALUGridElementType { cube, simplex };
  //! \brief available refinement types for ALUGrid
  enum ALUGridRefinementType { conforming, nonconforming };

  template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType,
      class Comm =
#if HAVE_MPI
        MPI_Comm
#else
        No_Comm
#endif
      >
  class ALUGrid;
}
#endif
