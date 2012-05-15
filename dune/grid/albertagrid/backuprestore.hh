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

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,path,fileprefix)  */
    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      std::string filename ( path );
      if( fileprefix != "" )
        filename += "/" + fileprefix ;

      grid.writeGridXdr( filename, 0.0 );
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)
        \note This method is not available for AlbertGrid.
              Use try/catch to catch the NotImplemented exception
              and fall back to the other backup method. */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore using streams not implemented." );
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(path,fileprefix) */
    static Grid *restore ( const std::string &path, const std::string &fileprefix )
    {
      std::string filename ( path );
      if( fileprefix != "" )
        filename += "/" + fileprefix ;

      Grid *grid = new Grid;
      double time; // ignore time
      grid->readGridXdr( filename, time );
      return grid;
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream)
        \note This method is not available for AlbertGrid.
              Use try/catch to catch the NotImplemented exception
              and fall back to the other restore method. */
    static Grid *restore ( std::istream &stream )
    {
      DUNE_THROW( NotImplemented, "backup / restore using streams not implemented." );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALBERTAGRID_BACKUPRESTORE_HH
