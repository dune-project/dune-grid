// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFGEOGRID_HH
#define DUNE_DGFGEOGRID_HH

#include <dune/grid/geogrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

namespace Dune
{

  /************************************************************************
  * Warning:
  * Reading DGF files directly into a GeometryGrid is a dirty hack for
  * two reasons:
  *   1) The host grid and coordinate function are never deleted (dangling
  *      pointers).
  *   2) The coordinate function has to provide a default constructor
  ************************************************************************/


  // MacroGrid :: Impl for GeomegryGrid
  // ----------------------------------

  /** \cond */
  template< class HostGrid, class CoordFunction >
  struct MacroGrid :: Impl< GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef MPIHelper :: MPICommunicator MPICommunicator;

    inline static GeometryGrid< HostGrid, CoordFunction > *
    generate ( MacroGrid &macroGrid,
               const char *filename,
               MPICommunicator communicator = MPIHelper :: getCommunicator() );
  };


  template< class HostGrid, class CoordFunction >
  inline GeometryGrid< HostGrid, CoordFunction > *
  MacroGrid :: Impl< GeometryGrid< HostGrid, CoordFunction > >
  :: generate ( MacroGrid &macroGrid,
                const char *filename,
                MPICommunicator communicator )
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;
    typedef MacroGrid :: Impl< HostGrid > HostImpl;

    HostGrid *hostGrid = HostImpl :: generate( macroGrid, filename, communicator );
    assert( hostGrid != 0 );
    const CoordFunction *coordFunction = new CoordFunction;
    return new Grid( *hostGrid, *coordFunction );
  }



  // DGFGridInfo for GeometryGrid
  // ----------------------------

  template< class HostGrid, class CoordFunction >
  struct DGFGridInfo< GeometryGrid< HostGrid, CoordFunction > >
  {
    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HostGrid > :: refineStepsForHalf();
    }

    static double refineWeight ()
    {
      return -1.0;
    }
  };
  /** \endcond */

}

#endif
