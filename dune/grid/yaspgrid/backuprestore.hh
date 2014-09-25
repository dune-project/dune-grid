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
  template<int dim, class ctype>
  struct BackupRestoreFacility<YaspGrid<dim,EquidistantCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,EquidistantCoordinates<ctype,dim> > Grid;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      std::ostringstream filename_str;
      filename_str << filename;
      std::ofstream file( filename_str.str() );
      if( file )
      {
        // only write something if this is rank 0.
        if (grid.comm().rank() == 0)
          backup(grid,file);
        file.close();
      }
      else
        std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename_str.str() << "'" << std::endl;
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      stream << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: " << (grid.getRefineOption() ? "1" : "0") << std::endl;
      grid.begin()->coords.print(stream);
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore ( const std::string &filename )
    {
      std::ostringstream filename_str;
      filename_str << filename;
      std::ifstream file(filename_str.str());
      if( file )
        return restore(file);
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename_str.str() << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore ( std::istream &stream )
    {
      std::string input;

      int refinement;
      stream >> input >> input;
      stream >> refinement;
      std::cout << "Refinement level: " << refinement << std::endl;

      std::bitset<dim> periodic;
      bool b;
      stream >> input;
      for (int i=0; i<dim; i++)
      {
        stream >> b;
        periodic[i] = b;
      }
      std::cout << "Periodicity: " << periodic << std::endl;

      int overlap;
      stream >> input;
      stream >> overlap;
      std::cout << "Overlap size: " << overlap << std::endl;

      bool physicalOverlapSize;
      stream >> input;
      stream >> physicalOverlapSize;
      std::cout << "Keep physical overlap size: " << physicalOverlapSize << std::endl;

      Dune::FieldVector<ctype,dim> h;
      stream >> input >> input >> input >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> h[i];
      std::cout << "Meshsize: " << h << std::endl;

      Dune::array<int,dim> s;
      stream >> input;
      for (int i=0; i<dim; i++)
        stream >> s[i];
      std::cout << "Size: " << s << std::endl;

      // the constructor takes the upper right corner...
      Dune::FieldVector<ctype,dim> length(h);
      for (int i=0; i<dim; i++)
        h[i] *= s[i];

#if HAVE_MPI
      Grid* grid = new Dune::YaspGrid<dim>(MPI_COMM_WORLD,length, s, periodic, overlap);
#else
      Grid* grid = new Dune::YaspGrid<dim>(length, s, periodic, overlap);
#endif

      grid->refineOptions(physicalOverlapSize);
      grid->globalRefine(refinement);

      return grid;
    }
  };

  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class ctype>
  struct BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,TensorProductCoordinates<ctype,dim> > Grid;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      std::ostringstream filename_str;
      filename_str << filename << grid.comm().rank();
      std::ofstream file( filename_str.str() );
      if( file )
      {
        backup(grid,file);
        file.close();
      }
      else
        std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename_str.str() << "'" << std::endl;
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      stream << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: " << (grid.getRefineOption() ? "1" : "0") << std::endl;
      grid.begin()->coords.print(stream);
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore ( const std::string &filename )
    {
      std::ostringstream filename_str;
      filename_str << filename;
      // TODO replace this by collective communication (handles sequential grids too)
#if HAVE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      filename_str << rank;
#endif
      std::ifstream file(filename_str.str());
      if( file )
        return restore(file);
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename_str.str() << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore ( std::istream &stream )
    {
      std::string input;

      int refinement;
      stream >> input >> input;
      stream >> refinement;
      std::cout << "Refinement level: " << refinement << std::endl;

      std::bitset<dim> periodic;
      bool b;
      stream >> input;
      for (int i=0; i<dim; i++)
      {
        stream >> b;
        periodic[i] = b;
      }
      std::cout << "Periodicity: " << periodic << std::endl;

      int overlap;
      stream >> input;
      stream >> overlap;
      std::cout << "Overlap size: " << overlap << std::endl;

      bool physicalOverlapSize;
      stream >> input;
      stream >> physicalOverlapSize;
      std::cout << "Keep physical overlap size: " << physicalOverlapSize << std::endl;

      Dune::array<int,dim> offset;
      stream >> input;
      for (int i=0; i<dim; i++)
        stream >> offset[i];
      std::cout << "ProcessorOffset: " << offset << std::endl;

      Dune::array<std::vector<ctype>,dim> coords;
      stream >> input >> input >> input >> input;
      for (int d=0; d<dim; d++)
      {
        stream >> input >> input;
        int size;
        stream >> size;
        stream >> input;
        ctype tmp;
        coords[d].resize(size);
        for (int i=0; i<size; i++)
        {
          stream >> tmp;
          coords[d][i] = tmp;
        }
      }

      // TODO this is a stub. It does treat the grid on the processor as a new grid.
      Grid* grid = new Grid(MPI_COMM_WORLD, coords, periodic, overlap);

      grid->refineOptions(physicalOverlapSize);
      grid->globalRefine(refinement);

      return grid;
    }
  };
} // namespace Dune

#endif // #ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
