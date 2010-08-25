// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CAPABILITIES_HH
#define DUNE_CAPABILITIES_HH

/** \file
    \brief A set of traits classes to store static information about grid implementation
 */

namespace Dune
{

  /** \brief Contains all capabilities classes */
  namespace Capabilities
  {

    /** \brief Specialize with 'true' for all codims that a grid implements entities for. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid, int codim>
    struct hasEntity
    {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if implementation supports parallelism. (default=false)
        \ingroup GICapabilities
     */
    template<class Grid>
    struct isParallel
    {
      static const bool v = false;
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

        \ingroup GICapabilities
     */
    template <class Grid>
    struct threadSafe {
      static const bool v = false;
    };

    /** \brief Specialize with 'true' if the grid implementation is thread safe, while it is not modified. (default=false)

        This capability is 'true' if the grid thread safe, <b>at least</b> while it is not modified.
        This means especially that all grid modification, like load balancing and
        grid adaption are not necessarily thread safe.

        \sa threadSafe

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

    template<class Grid, int codim>
    struct hasEntity<const Grid, codim>
    {
      static const bool v = Dune::Capabilities::hasEntity<Grid,codim>::v;
    };

    template<class Grid>
    struct isParallel<const Grid>
    {
      static const bool v = Dune::Capabilities::isParallel<Grid>::v;
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

#endif // DUNE_CAPABILITIES_HH
