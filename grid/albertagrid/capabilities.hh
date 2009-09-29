// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_CAPABILITIES_HH
#define DUNE_ALBERTA_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    /** \brief   AlbertaGrid has entities for all codimensions
     *  \ingroup AlbertaGrid
     */
    template< int dim, int dimworld, int codim >
    struct hasEntity< AlbertaGrid< dim, dimworld >, codim >
    {
      static const bool v = true;
    };

    /** \brief   AlbertaGrid is not levelwise conforming
     *           (since it uses bisection)
     *  \ingroup AlbertaGrid
     */
    template< int dim, int dimworld >
    struct isLevelwiseConforming< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    /** \brief   AlbertaGrid is leafwise conforming
     * \ingroup AlbertaGrid
     */
    template< int dim, int dimworld >
    struct isLeafwiseConforming< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    /** \brief AlbertaGrid does not support hanging nodes
     *  \ingroup AlbertaGrid
     */
    template< int dim, int dimworld >
    struct hasHangingNodes< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    /** \brief   AlbertaGrid has backup and restore facilities
     *  \ingroup AlbertaGrid
     */
    template< int dim, int dimworld >
    struct hasBackupRestoreFacilities< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };



    // non-standard capabilities
    // -------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< int dim, int dimworld >
    struct hasHierarchicIndexSet< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

  }

}

#endif // #ifndef DUNE_ALBERTA_CAPABILITIES_HH
