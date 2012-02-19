// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_CAPABILITIES_HH
#define DUNE_ALUGRID_CAPABILITIES_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/common/checkparallel.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

/** @file
 *  @author Robert Kloefkorn
 *  @brief Capabilities for ALUGrid
 */

namespace Dune
{

  namespace Capabilities
  {

    // Capabilities for ALUGrid
    // -------------------------------

    /** \brief ALUGrid has only one geometry type for codim 0 entities
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct hasSingleGeometryType< ALUGrid< dim, dimworld, eltype, refinementtype, Comm > >
    {
      static const bool v = true;
      static const unsigned int topologyId = (eltype == cube) ?
                                             GenericGeometry :: CubeTopology< dim > :: type :: id :
                                             GenericGeometry :: SimplexTopology< dim > :: type :: id ;
    };

    /** \brief ALUGrid has entities for all codimension
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm, int cdim >
    struct hasEntity< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUGrid is parallel when Comm == MPI_Comm
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype >
    struct isParallel< ALUGrid< dim, dimworld, eltype, refinementtype, No_Comm > >
    {
      static const bool v = false;
    };

    /** \brief ALUGrid is parallel when Comm == MPI_Comm
       \ingroup ALUGrid
     */
#if ALU3DGRID_PARALLEL
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype >
    struct isParallel< ALUGrid< 3, 3, eltype, refinementtype,  MPI_Comm > >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUGrid can communicate when Comm == MPI_Comm
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, int codim >
    struct canCommunicate< ALUGrid< dim, dimworld, eltype, refinementtype, No_Comm >, codim >
    {
      static const bool v = false;
    };

    /** \brief ALUGrid can communicate
       \ingroup ALUGrid
     */
#if ALU3DGRID_PARALLEL
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, int codim >
    struct canCommunicate< ALUGrid< 3, 3, eltype, refinementtype, MPI_Comm >, codim >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUGrid has conforming level grids
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct isLevelwiseConforming< ALUGrid< dim, dimworld, eltype, refinementtype, Comm > >
    {
      static const bool v = refinementtype == nonconforming;
    };

    /** \brief ALUGrid has conforming level grids
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct isLeafwiseConforming< ALUGrid< dim, dimworld, eltype, refinementtype, Comm > >
    {
      static const bool v = refinementtype == conforming ;
    };

    /** \brief ALUGrid has backup and restore facilities
       \ingroup ALUGrid
     */
    template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct hasBackupRestoreFacilities< ALUGrid< dim, dimworld, eltype, refinementtype, Comm > >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

} //end  namespace Dune

#endif // #ifdef ENABLE_ALUGRID

#endif // #ifdef DUNE_ALUGRID_CAPABILITIES_HH
