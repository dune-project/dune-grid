// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_CAPABILITIES_HH
#define DUNE_ALU2DGRID_CAPABILITIES_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

#include <dune/grid/alugrid/2d/alu2dinclude.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>


/** @file
 *  @author Robert Kloefkorn
 *  @brief Capabilities for 2d ALUGrids
 */

namespace Dune
{

  template< int dim, int dimworld >
  class ALUSimplexGrid;

  template< int dim, int dimw >
  class ALUCubeGrid;

  template< int dim, int dimworld >
  class ALUConformGrid;



  namespace Capabilities
  {


    // Capabilities for ALUSimplexGrid
    // -------------------------------

    /** \brief ALUSimplexGrid has only one geometry type for codim 0 entities
       \ingroup ALUSimplexGrid
     */
    template< int dimworld >
    struct hasSingleGeometryType< ALUSimplexGrid< 2, dimworld > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: SimplexTopology< 2 > :: type :: id ;
    };


    /** \brief ALUSimplexGrid has entities for all codimension
       \ingroup ALUSimplexGrid
     */
    template< int dimworld, int cdim >
    struct hasEntity< ALUSimplexGrid< 2, dimworld >, cdim >
    {
      static const bool v = true;
    };

#if ALU2DGRID_PARALLEL
    /** \brief ALUSimplexGrid is parallel
       \ingroup ALUSimplexGrid
     */
    //- default is false
    template< int dimworld >
    struct isParallel< ALUSimplexGrid< 2, dimworld > >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

#if ALU2DGRID_PARALLEL
    /** \brief ALUSimplexGrid can communicate
       \ingroup ALUSimplexGrid
     */
    //- default is false
    template< int dimworld >
    struct canCommunicate< ALUSimplexGrid< 2, dimworld >, 0 >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

    /** \brief ALUSimplexGrid has conforming level grids
       \ingroup ALUSimplexGrid
     */
    template< int dimworld >
    struct isLevelwiseConforming< ALUSimplexGrid< 2, dimworld > >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid has backup and restore facilities
       \ingroup ALUSimplexGrid
     */
    template< int dimworld >
    struct hasBackupRestoreFacilities< ALUSimplexGrid< 2, dimworld > >
    {
      static const bool v = true;
    };



    // Capabilities for ALUCubeGrid
    // ----------------------------

    /** \brief ALUCubeGrid has only one geometry type for codim 0 entities
       \ingroup ALUCubeGrid
     */
    template< int wdim >
    struct hasSingleGeometryType< ALUCubeGrid< 2, wdim > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< 2 > :: type :: id ;
    };

    /** \brief ALUCubeGrid has entities for all codimension
       \ingroup ALUCubeGrid
     */
    template< int wdim, int cdim >
    struct hasEntity< Dune::ALUCubeGrid< 2, wdim >, cdim >
    {
      static const bool v = true;
    };

#if ALU2DGRID_PARALLEL
    /** \brief ALUCubeGrid is parallel
       \ingroup ALUCubeGrid
     */
    //- default is false
    template< int dimworld >
    struct isParallel< ALUCubeGrid< 2, dimworld > >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

#if ALU2DGRID_PARALLEL
    /** \brief ALUCubeGrid can communicate
       \ingroup ALUCubeGrid
     */
    //- default is false
    template< int dimworld >
    struct canCommunicate< ALUCubeGrid< 2, dimworld >, 0 >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

    /** \brief ALUCubeGrid has conforming level grids
       \ingroup ALUCubeGrid
     */
    template<int wdim>
    struct isLevelwiseConforming< Dune::ALUCubeGrid< 2, wdim > >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid has backup and restore facilities
       \ingroup ALUCubeGrid
     */
    template<int wdim>
    struct hasBackupRestoreFacilities< Dune::ALUCubeGrid< 2, wdim > >
    {
      static const bool v = true;
    };



    // Capabilities for ALUConformGrid
    // -------------------------------

    /** \brief ALUConformGrid has only one geometry type for codim 0 entities
       \ingroup ALUConformGrid
     */
    template< int dimworld >
    struct hasSingleGeometryType< ALUConformGrid< 2, dimworld > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: SimplexTopology< 2 > :: type :: id ;
    };

    /** \brief ALUConformGrid has entities for all codimension
       \ingroup ALUConformGrid
     */
    template< int dimworld, int cdim >
    struct hasEntity< ALUConformGrid< 2, dimworld >, cdim >
    {
      static const bool v = true;
    };

#if ALU2DGRID_PARALLEL
    /** \brief ALUConformGrid is parallel
       \ingroup ALUConformGrid
     */
    //- default is false
    template< int dimworld >
    struct isParallel< ALUConformGrid< 2, dimworld > >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

#if ALU2DGRID_PARALLEL
    /** \brief ALUConformGrid can communicate
       \ingroup ALUConformGrid
     */
    //- default is false
    template< int dimworld >
    struct canCommunicate< ALUConformGrid< 2, dimworld >, 0 >
    {
      static const bool v = true;
    };
#endif // #if ALU2DGRID_PARALLEL

    /** \brief ALUConformGrid has a conforming leaf grid
       \ingroup ALUConformGrid
     */
    template< int dimworld >
    struct isLeafwiseConforming< ALUConformGrid< 2, dimworld > >
    {
      static const bool v = true;
    };

    /** \brief ALUConformGrid has backup and restore facilities
       \ingroup ALUConformGrid
     */
    template< int dimworld >
    struct hasBackupRestoreFacilities< ALUConformGrid< 2, dimworld > >
    {
      static const bool v = true;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifdef ENABLE_ALUGRID

#endif // #ifndef DUNE_ALU2DGRID_CAPABILITIES_HH
