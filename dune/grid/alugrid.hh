// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// only include this code, if HAVE_ALUGRID is true
#if HAVE_ALUGRID

#include <dune/grid/alugrid/3d/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

// 2d version
#include <dune/grid/alugrid/2d/alugrid.hh>
#include <dune/grid/alugrid/2d/alu2dgridfactory.hh>

#include <dune/grid/alugrid/common/persistentcontainer.hh>

#include <alugrid_defineparallel.h>
#if ALU3DGRID_BUILD_FOR_PARALLEL && ! HAVE_MPI
#error When using parallel ALUGrid, you need to
#error * configure dune with --enable-parallel
#error * add the DUNEMPI*FLAGS
#endif

/** @file
    @author Robert Kloefkorn
    @brief Provides base classes for ALUGrid
 **/

namespace Dune
{

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 3D grid with support for hexahedrons.
     @ingroup GridImplementations
     @ingroup ALUCubeGrid
     The ALUCubeGrid implements the Dune GridInterface for 3d hexahedral meshes.
     This grid can be locally adapted (non-conforming) and used in parallel
     computations using dynamic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports hexahedrons - a 2d/3d simplex
     grid is also available via the grid implementation ALUSimplexGrid or ALUConformGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)

     \li Available Implementations
          - Dune::ALUCubeGrid<3,3>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template< int dim, int dimworld >
  class ALUCubeGrid;



  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUSimplexGrid
     The ALUSimplexGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted (non-conforming) and used in parallel
     computations using dynamic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports triangular/tetrahedral elements -
     a 3d cube
     grid is also available via the grid implementation ALUCubeGrid or ALUConformGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)

     \li Available Implementations
            - Dune::ALUSimplexGrid<3,3>
            - Dune::ALUSimplexGrid<2,2>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template< int dim, int dimworld >
  class ALUSimplexGrid;


  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUConformGrid
     The ALUConformGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted (conforming) and used in parallel
     computations using dynamic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports triangular/tetrahedral elements -
     a 3d cube
     grid is also available via the grid implementation ALUCubeGrid or ALUSimplexGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     \li Available Implementations
            - Dune::ALUConformGrid<2,2>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template <int dim, int dimworld>
  class ALUConformGrid;

} //end  namespace Dune

#endif // #ifdef HAVE_ALUGRID

#endif
