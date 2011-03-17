// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ALBERTAGRID_BACKUPRESTORE_HH
#define DUNE_GRID_ALBERTAGRID_BACKUPRESTORE_HH

#include <dune/grid/common/backuprestore.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int >
  class AlbertaGrid;



  // BackupRestoreFacility for AlbertaGrid
  // -------------------------------------

  template< int dim, int dimworld >
  struct BackupRestoreFacility< AlbertaGrid< dim, dimworld > >
  {
    typedef AlbertaGrid< dim, dimworld > Grid;

    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      const std::string filename( path + "/" + fileprefix );
      return grid.writeXdr( filename, 0.0 );
    }

    static void backup ( const Grid &grid, const std::ostream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore using streams not implemented." );
    }

    static Grid *restore ( const std::string &path, const std::string &fileprefix )
    {
      const std::string filename( path + "/" + fileprefix );
      Grid *grid = new Grid;
      double time; // ignore time
      grid->readGridXdr( filename, time );
      return grid;
    }

    static Grid *restore ( const std::istream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore using streams not implemented." );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALBERTAGRID_BACKUPRESTORE_HH
