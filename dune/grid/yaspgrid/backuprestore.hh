// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
#define DUNE_GRID_YASPGRID_BACKUPRESTORE_HH

//- system headers
#include <fstream>

//- Dune headers
#include <dune/common/exceptions.hh>
#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/alugrid/common/declaration.hh>

namespace Dune
{

  /** \copydoc Dune::BackupRestoreFacility */
  template<class ctype, int dim>
  struct BackupRestoreFacility<YaspGrid<dim,EquidistantCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,EquidistantCoordinates<ctype,dim> > Grid;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      std::ofstream file( filename.c_str() );
      if( file )
      {
        backup(grid,file);
      }
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename << "'" << std::endl;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      stream << "Refinement level: " << grid.maxLevel() << std::endl;
      grid.begin()->coords.print(stream);

    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore ( const std::string &filename )
    {
      std::ifstream file( filename.c_str() );
      if( file )
      {
        //TODO
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
      // TODO
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
