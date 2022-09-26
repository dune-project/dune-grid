// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_CAPABILITIES_HH
#define DUNE_GRID_COMMON_CAPABILITIES_HH

/** \file
    \brief A set of traits classes to store static information about grid implementation
 */

namespace Dune
{

  /** \brief Contains all capabilities classes */
  namespace Capabilities
  {

    /** \brief Specialize with 'true' for if the codimension 0 entity
        of the grid has only one possible geometry type. In this case the
        topologyId of this geometry type has also to be specified.
        (default=false, topologyId=undefined)
        \ingroup GICapabilities
     */
    template<class Grid>
    struct hasSingleGeometryType
    {
      static const bool v = false;
      // this value will be initialized with something big
      // since it is invalid
      static const unsigned int topologyId = ~0u;
    };

    /** \brief Specialize with 'true' if the grid is a Cartesian grid.
        Cartesian grids satisfy the following properties:
          - all geometries are axis-aligned hypercubes
          - The unit outer normal for the i-th intersection
            can be computed by the following code:
          \code
             FieldVector< ctype, dim > n( 0 );
             n[ i / 2 ] = ctype( 2*(i % 2) - 1 );
          \endcode
        (default=false).
        \ingroup GICapabilities
     */
    template<class Grid>
    struct isCartesian
    {
      // default value is false
      static const bool v = false;
    };

    /** \brief Specialize with 'true' for all codims that a grid implements entities for. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid, int codim>
    struct hasEntity
    {
      static const bool v = false;
    };

    /**
     * \brief specialize with 'true' for all codims that a grid provides an iterator for (default=hasEntity<codim>::v)
     *
     * \note Being able to iterate over a codimension implies that the grid
     *       provides entities for that codimension.
     * \note Any fully conforming DUNE grid implementation must be able to
     *       iterate over all codim 0 entities (i.e., elements).
     *
     * \ingroup GICapabilities
     **/
    template< class Grid, int codim >
    struct hasEntityIterator
    {
      static const bool v = hasEntity<Grid, codim>::v;
    };

    /** \brief Specialize with 'false' for all codims that a grid does not
      implement geometries for. (default=true)
        \ingroup GICapabilities
     */
    template<class Grid, int codim>
    struct hasGeometry
    {
      static const bool v = true;
    };

    /** \brief specialize with 'true' for all codims that a grid can communicate data on (default=false)
     *
     *  \note Being able to communicate data on a codimension implies that the
     *        grid provides entities for that codimension.
     *
     *  \ingroup GICapabilities
     */
    template< class Grid, int codim >
    struct canCommunicate
    {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if implementation guarantees conforming level grids. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid>
    struct isLevelwiseConforming
    {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if implementation guarantees a conforming leaf grid. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid>
    struct isLeafwiseConforming
    {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if implementation provides backup and restore facilities. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid>
    struct hasBackupRestoreFacilities
    {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if the grid implementation is thread safe. (default=false)

        This capability is 'true' if the grid is <b>always</b> thread safe.
        This means especially that all grid modification, like load balancing and
        grid adaption are also thread safe.

        \sa viewThreadSafe

        \note that the communicate method can only be called by one individual thread,
        as the whole Dune parallel components are not (i.e. cannot be) thread safe.

        \ingroup GICapabilities
     */
    template <class Grid>
    struct threadSafe {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if the grid implementation is thread safe, while it is not modified. (default=false)

        This capability is 'true' if the grid is thread safe, <b>at least</b> while working on a GridView.
        This means especially that all grid modification, like load balancing and
        grid adaption are not necessarily thread safe.

        \sa threadSafe

        \note that the communicate method can only be called by one individual thread,
        as the whole Dune parallel components are (i.e. cannot be) not thread safe.

        \note the methods leafGridView(), levelGridView(level) on the Grid can only be called single-threaded

        \note calling the methods indexSet(), idSet(), globalIdSet() on the Grid or the GridView is only allowed,
        if they were called once before starting the threads.

        \code

        \endcode

        \ingroup GICapabilities
     */
    template <class Grid>
    struct viewThreadSafe {
      static const bool v = false;
    };

    /*
       forward
       Capabilities::Something<const Grid>
       to
       Capabilities::Something<Grid>
     */

    template<class Grid>
    struct hasSingleGeometryType< const Grid >
    {
      static const bool v = Dune::Capabilities::hasSingleGeometryType<Grid>::v;
      static const unsigned int topologyId =
        Dune::Capabilities::hasSingleGeometryType<Grid>::topologyId;
    };

    template<class Grid>
    struct isCartesian< const Grid >
    {
      static const bool v = Dune::Capabilities::isCartesian<Grid>::v;
    };

    template<class Grid, int codim>
    struct hasEntity<const Grid, codim>
    {
      static const bool v = Dune::Capabilities::hasEntity<Grid,codim>::v;
    };

    template< class Grid, int codim >
    struct hasEntityIterator< const Grid, codim >
    {
      static const bool v = Dune::Capabilities::hasEntityIterator< Grid, codim >::v;
    };

    template< class Grid, int codim >
    struct canCommunicate< const Grid, codim >
    {
      static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
    };

    template<class Grid>
    struct isLevelwiseConforming<const Grid>
    {
      static const bool v = Dune::Capabilities::isLevelwiseConforming<Grid>::v;
    };

    template<class Grid>
    struct isLeafwiseConforming<const Grid>
    {
      static const bool v = Dune::Capabilities::isLeafwiseConforming<Grid>::v;
    };

    template<class Grid>
    struct hasBackupRestoreFacilities<const Grid>
    {
      static const bool v = Dune::Capabilities::hasBackupRestoreFacilities<Grid>::v;
    };

    template <class Grid>
    struct threadSafe<const Grid> {
      static const bool v = Dune::Capabilities::threadSafe<Grid>::v;
    };

    template <class Grid>
    struct viewThreadSafe<const Grid> {
      static const bool v = Dune::Capabilities::viewThreadSafe<Grid>::v;
    };

  }

}

#endif // DUNE_GRID_COMMON_CAPABILITIES_HH
