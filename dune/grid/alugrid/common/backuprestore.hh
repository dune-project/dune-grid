// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
#define DUNE_GRID_ALUGRID_BACKUPRESTORE_HH

//- system headers
#include <fstream>

//- Dune headers
#include <dune/common/exceptions.hh>
#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/alugrid/common/declaration.hh>

namespace Dune
{

  /** \copydoc Dune::BackupRestoreFacility */
  template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
  struct BackupRestoreFacility< ALUGrid< dim, dimworld, elType, refineType, Comm > >
  {
    // type of grid
    typedef ALUGrid< dim, dimworld, elType, refineType, Comm > Grid;

    static std::string createFilename( const std::string &path, const std::string &fileprefix )
    {
      std::string filename( path );
      if( fileprefix.size() > 0 )
      {
        filename += "/" + fileprefix ;
      }
      else if( filename[ filename.size() - 1 ] == char('/') )
      {
        filename += "/alugrid";
      }
      return filename;
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,path,fileprefix)  */
    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      std::string filename ( createFilename( path, fileprefix ) );
      std::ofstream file( filename.c_str() );
      if( file )
      {
        // call backup on grid
        backup( grid, file );
      }
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename << "'" << std::endl;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      // call backup on grid
      grid.backup( stream );
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(path,fileprefix) */
    static Grid *restore ( const std::string &path, const std::string &fileprefix )
    {
      // Problem here: how to pass boundary projections
      std::string filename( createFilename( path, fileprefix ) );
      std::ifstream file( filename.c_str() );
      if( file )
      {
        return restore( file );
      }
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore ( std::istream &stream )
    {
      // Problem here: how to pass boundary projections
      Grid* grid = new Grid();
      grid->restore( stream );
      return grid;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
