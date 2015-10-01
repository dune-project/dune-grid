// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
#define DUNE_GEOGRID_BACKUPRESTORE_HH

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

    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      // notice: We should also backup the coordinate function
      HostBackupRestoreFacility::backup( grid.hostGrid(), path, fileprefix );
    }

    static void backup ( const Grid &grid, const std::ostream &stream )
    {
      // notice: We should also backup the coordinate function
      HostBackupRestoreFacility::backup( grid.hostGrid(), stream );
    }

    static Grid *restore ( const std::string &path, const std::string &fileprefix )
    {
      // notice: We should also restore the coordinate function
      HostGrid *hostGrid = HostBackupRestoreFacility::restore( path, fileprefix );
      CoordFunction *coordFunction = new CoordFunction();
      return new Grid( hostGrid, coordFunction );
    }

    static Grid *restore ( const std::istream &stream )
    {
      // notice: We should also restore the coordinate function
      HostGrid *hostGrid = HostBackupRestoreFacility::restore( stream );
      CoordFunction *coordFunction = new CoordFunction();
      return new Grid( hostGrid, coordFunction );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
