// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
#define DUNE_GEOGRID_BACKUPRESTORE_HH

#include <dune/grid/utility/grapedataioformattypes.hh>

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

    public:
      template< GrapeIOFileFormatType type >
      bool writeGrid ( const std::string &filename, double time ) const
      {
        return asImp().hostGrid().template writeGrid< type >( filename, time );
      }

      template< GrapeIOFileFormatType type >
      bool readGrid ( const std::string &filename, double &time )
      {
        const bool success
          = asImp().hostGrid().template readGrid< type >( filename, time );
        asImp().update();
        return success;
      }

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

  }

}

#endif // #ifndef DUNE_GEOGRID_BACKUPRESTORE_HH
