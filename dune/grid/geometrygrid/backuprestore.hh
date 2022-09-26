// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
#define DUNE_GEOGRID_BACKUPRESTORE_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/backuprestore.hh>

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/capabilities.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // BackupRestoreFacilities
    // -----------------------

    template< class Grid, bool hasBackupRestoreFacilities = Capabilities::hasBackupRestoreFacilities< Grid > ::v >
    class BackupRestoreFacilities
    {};

    template< class Grid >
    class BackupRestoreFacilities< Grid, true >
    {
      typedef BackupRestoreFacilities< Grid, true > This;

    protected:
      BackupRestoreFacilities ()
      {}

    private:
      BackupRestoreFacilities ( const This & );
      This &operator= ( const This & );

    protected:
      const Grid &asImp () const
      {
        return static_cast< const Grid & >( *this );
      }

      Grid &asImp ()
      {
        return static_cast< Grid & >( *this );
      }
    };

  } // namespace GeoGrid



  // BackupRestoreFacility for GeometryGrid
  // --------------------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  struct BackupRestoreFacility< GeometryGrid< HostGrid, CoordFunction, Allocator > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;
    typedef BackupRestoreFacility< HostGrid > HostBackupRestoreFacility;

    /// \brief Backup the grid to file or stream
    template <class Output>
    static void backup ( const Grid &grid, const Output &filename_or_stream )
    {
      // notice: We should also backup the coordinate function
      HostBackupRestoreFacility::backup( grid.hostGrid(), filename_or_stream );
    }

    /// \brief Restore the grid from file or stream
    template <class Input>
    static Grid *restore ( const Input &filename_or_stream )
    {
      // notice: assumes the CoordFunction to be default-constructible
      return restore_impl(filename_or_stream, std::is_default_constructible<CoordFunction>{});
    }

  private:
    template <class Input>
    static Grid *restore_impl ( const Input &filename_or_stream, std::true_type )
    {
      // notice: We should also restore the coordinate function
      HostGrid *hostGrid = HostBackupRestoreFacility::restore( filename_or_stream );
      CoordFunction *coordFunction = new CoordFunction();
      return new Grid( hostGrid, coordFunction );
    }

    template <class Input>
    static Grid *restore_impl ( const Input &filename_stream, std::false_type )
    {
      DUNE_THROW(NotImplemented,
        "Restoring a GeometryGrid with a CoordFunction that is not default-constructible is not implemented.");
      return nullptr;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
