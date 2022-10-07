// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_CAPABILITIES_HH
#define DUNE_ALBERTA_CAPABILITIES_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_ALBERTA

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

    /** \brief AlbertaGrid has only one geometry type for codim 0 entities
       \ingroup AlbertaGrid
     */
    template< int dim, int dimworld >
    struct hasSingleGeometryType< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::simplex(dim).id();
    };


    /** \brief   AlbertaGrid has entities for all codimensions
     *  \ingroup AlbertaGrid
     */
    template< int dim, int dimworld, int codim >
    struct hasEntity< AlbertaGrid< dim, dimworld >, codim >
    {
      static const bool v = true;
    };

    /**
     * \brief AlbertaGrid can iterate over all codimensions
     *
     * \ingroup AlbertaGrid
     **/
    template< int dim, int dimworld, int codim >
    struct hasEntityIterator< AlbertaGrid< dim, dimworld >, codim >
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

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_CAPABILITIES_HH
