// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_CAPABILITIES_HH
#define DUNE_ALU3DGRID_CAPABILITIES_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/3d/alu3dinclude.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>


/** @file
 *  @author Robert Kloefkorn
 *  @brief Capabilities for 3d ALUGrids
 */

namespace Dune
{

  template< int dim, int dimworld >
  class ALUCubeGrid;

  template< int dim, int dimworld >
  class ALUSimplexGrid;


  namespace Capabilities
  {

    // Capabilities for ALUCubeGrid
    // ----------------------------

    /** \struct isLeafwiseConforming
       \ingroup ALUCubeGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUCubeGrid
     */

    /** \brief ALUCubeGrid has only one geometry type for codim 0 entities
       \ingroup ALUCubeGrid
     */
    template< >
    struct hasSingleGeometryType< ALUCubeGrid< 3, 3 > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< 3 > :: type :: id ;
    };


    /** \brief ALUCubeGrid has entities for all codimension
       \ingroup ALUCubeGrid
     */
    template< int cdim >
    struct hasEntity< ALUCubeGrid< 3, 3 >, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid is parallel
       \ingroup ALUCubeGrid
     */
#if ALU3DGRID_PARALLEL
    template<>
    struct isParallel< ALUCubeGrid< 3, 3 > >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUCubeGrid can communicate
       \ingroup ALUCubeGrid
     */
#if ALU3DGRID_PARALLEL
    template< int codim >
    struct canCommunicate< ALUCubeGrid< 3, 3 >, codim >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUCubeGrid has conforming level grids
       \ingroup ALUCubeGrid
     */
    template<>
    struct isLevelwiseConforming< ALUCubeGrid< 3, 3 > >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid has backup and restore facilities
       \ingroup ALUCubeGrid
     */
    template<>
    struct hasBackupRestoreFacilities< ALUCubeGrid< 3, 3 > >
    {
      static const bool v = true;
    };



    // Capabilities for ALUSimplexGrid
    // -------------------------------

    /** \struct isLeafwiseConforming
       \ingroup ALUSimplexGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUSimplexGrid
     */

    /** \brief ALUSimplexGrid has only one geometry type for codim 0 entities
       \ingroup ALUSimplexGrid
     */
    template< >
    struct hasSingleGeometryType< ALUSimplexGrid< 3, 3 > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: SimplexTopology< 3 > :: type :: id ;
    };

    /** \brief ALUSimplexGrid has entities for all codimension
       \ingroup ALUSimplexGrid
     */
    template< int cdim >
    struct hasEntity< ALUSimplexGrid< 3, 3 >, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid is parallel
       \ingroup ALUSimplexGrid
     */
#if ALU3DGRID_PARALLEL
    template<>
    struct isParallel< ALUSimplexGrid< 3, 3 > >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUSimplexGrid can communicate
       \ingroup ALUSimplexGrid
     */
#if ALU3DGRID_PARALLEL
    template< int codim >
    struct canCommunicate< ALUSimplexGrid< 3, 3 >, codim >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUSimplexGrid has conforming level grids
       \ingroup ALUSimplexGrid
     */
    template<>
    struct isLevelwiseConforming< ALUSimplexGrid< 3, 3 > >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid has backup and restore facilities
       \ingroup ALUSimplexGrid
     */
    template<>
    struct hasBackupRestoreFacilities< ALUSimplexGrid< 3, 3 > >
    {
      static const bool v = true;
    };

    // Capabilities for ALUGrid
    // -------------------------------

    /** \struct isLeafwiseConforming
       \ingroup ALUGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUGrid
     */

    /** \brief ALUGrid has only one geometry type for codim 0 entities
       \ingroup ALUGrid
     */
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct hasSingleGeometryType< ALUGrid< 3, 3, eltype, refinementtype, Comm > >
    {
      static const bool v = true;
      static const unsigned int topologyId = (eltype == cube) ?
                                             GenericGeometry :: CubeTopology< 3 > :: type :: id :
                                             GenericGeometry :: SimplexTopology< 3 > :: type :: id ;
    };

    /** \brief ALUGrid has entities for all codimension
       \ingroup ALUGrid
     */
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm, int cdim >
    struct hasEntity< ALUGrid< 3, 3, eltype, refinementtype, Comm >, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUGrid is parallel
       \ingroup ALUGrid
     */
#if ALU3DGRID_PARALLEL
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct isParallel< ALUGrid< 3, 3, eltype, refinementtype, Comm > >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUGrid can communicate
       \ingroup ALUGrid
     */
#if ALU3DGRID_PARALLEL
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, int codim >
    struct canCommunicate< ALUGrid< 3, 3, eltype, refinementtype >, codim >
    {
      static const bool v = true;
    };
#endif

    /** \brief ALUGrid has conforming level grids
       \ingroup ALUGrid
     */
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct isLevelwiseConforming< ALUGrid< 3, 3, eltype, refinementtype, Comm > >
    {
      static const bool v = refinementtype == nonconforming;
    };

    /** \brief ALUGrid has conforming level grids
       \ingroup ALUGrid
     */
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype >
    struct isLeafwiseConforming< ALUGrid< 3, 3, eltype, refinementtype > >
    {
      static const bool v = refinementtype == conforming ;
    };

    /** \brief ALUGrid has backup and restore facilities
       \ingroup ALUGrid
     */
    template< ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
    struct hasBackupRestoreFacilities< ALUGrid< 3, 3, eltype, refinementtype, Comm > >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

} //end  namespace Dune

#endif // #ifdef ENABLE_ALUGRID

#endif
