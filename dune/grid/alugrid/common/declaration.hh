// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_FORWARDDECLARATION
#define DUNE_ALUGRID_FORWARDDECLARATION

namespace Dune {

  //! \brief basic element types for ALUGrid
  enum ALUGridElementType { cube, simplex };
  //! \brief available refinement types for ALUGrid
  enum ALUGridRefinementType { conforming, nonconforming };

  template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType >
  class ALUGrid;
}
#endif
