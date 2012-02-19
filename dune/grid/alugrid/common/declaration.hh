// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_FORWARDDECLARATION
#define DUNE_ALUGRID_FORWARDDECLARATION

//- common includes
#include <dune/common/collectivecommunication.hh>
#if HAVE_MPI
#include <dune/common/mpicollectivecommunication.hh>
#endif

//- ALUGrid includes
#include <dune/grid/alugrid/common/checkparallel.hh>

namespace Dune {

  //! \brief basic element types for ALUGrid
  enum ALUGridElementType { cube, simplex };
  //! \brief available refinement types for ALUGrid
  enum ALUGridRefinementType { conforming, nonconforming };

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for quadrilateral and hexahedral grid (template parameter cube)
     and simplicial meshes (template parameter simplex) in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUGrid

     The ALUGrid implements the Dune GridInterface for 2d quadrilateral and 3d hexahedral
     as well as 2d triangular and  3d tetrahedral meshes.
     This grid can be locally adapted (non-conforming and conforming bisection)
     and used in parallel computations using dynamic load balancing.

     @note
     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     \li Available Implementations
          - quadrilateral and hexahedral elements only nonconforming refinement
            - Dune::ALUGrid< 2, 2, cube, nonconforming >
            - Dune::ALUGrid< 2, 3, cube, nonconforming >
            - Dune::ALUGrid< 3, 3, cube, nonconforming >
          - simplicial elements and nonconforming refinement
            - Dune::ALUGrid< 2, 2, simplex, nonconforming >
            - Dune::ALUGrid< 2, 3, simplex, nonconforming >
            - Dune::ALUGrid< 3, 3, simplex, nonconforming >
          - simplicial elements and bisection refinement
            - Dune::ALUGrid< 2, 2, simplex, conforming >
            - Dune::ALUGrid< 2, 3, simplex, conforming >
            - Dune::ALUGrid< 3, 3, simplex, conforming > (work in progress)

     \note template parameter Comm defaults to MPI_Comm, if MPI is available, No_Comm  otherwise.

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType,
      class Comm =
#if ALU3DGRID_PARALLEL
        MPI_Comm
#else
        No_Comm
#endif
      >
  class ALUGrid;

  //- traits class for declaring base class for ALUGrid
  template <int dim, int dimw, ALUGridElementType elType, class Comm >
  struct ALUGridBaseGrid ;
}
#endif
