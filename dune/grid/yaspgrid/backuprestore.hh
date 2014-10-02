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
  /** \brief Implement load balancer that gets the info from a file
   * To backup and restore the load balancing information, the BackupRestoreFacility
   * needs a special load balancing object. Users dont need to touch this.
   */
  template<int d>
  class YLoadBalanceBackup : public YLoadBalance<d>
  {
  public:
    YLoadBalanceBackup(const Dune::array<int,d>& dims) : _dims(dims) {}

    virtual ~YLoadBalanceBackup() {}

    virtual void loadbalance(const Dune::array<int,d>&, int P, Dune::array<int,d>& dims) const
    {
      int prod = 1;
      for (int i=0; i<d; i++)
        prod *= _dims[i];
      if (P != prod)
        DUNE_THROW(Dune::Exception,"Your processor number doesn't match the grid to restore");
      dims = _dims;
    }

  private:
    Dune::array<int,d> _dims;
  };

  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class ctype>
  struct BackupRestoreFacility<YaspGrid<dim,EquidistantCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,EquidistantCoordinates<ctype,dim> > Grid;
    typedef typename Grid::Traits::CollectiveCommunication Comm;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      if (grid.comm().rank() == 0)
      {
        std::ofstream file(filename);
        if( file )
        {
          backup(grid,file);
          file.close();
        }
        else
          std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename << "'" << std::endl;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      stream << "Torus structure: ";
      for (int i=0; i<dim; i++)
        stream << grid.torus().dims(i) << " ";
      stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: " << (grid.getRefineOption() ? "1" : "0") << std::endl;
      stream << "Coarse Size: ";
      for (int i=0; i<dim; i++)
        stream << grid.levelSize(0,i) << " ";
      stream << std::endl;
      stream << "Meshsize: " ;
      for (int i=0; i<dim; i++)
        stream << grid.begin()->coords.meshsize(i,0) << " ";
      stream << std::endl;
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore (const std::string &filename, Comm comm = Comm())
    {
      std::ifstream file(filename);
      if( file )
        return restore(file,comm);
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore (std::istream &stream, Comm comm = Comm())
    {
      std::string input;

      Dune::array<int,dim> torus_dims;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> torus_dims[i];

      int refinement;
      stream >> input >> input;
      stream >> refinement;

      std::bitset<dim> periodic;
      bool b;
      stream >> input;
      for (int i=0; i<dim; i++)
      {
        stream >> b;
        periodic[i] = b;
      }

      int overlap;
      stream >> input;
      stream >> overlap;

      bool physicalOverlapSize;
      stream >> input;
      stream >> physicalOverlapSize;

      Dune::array<int,dim> coarseSize;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> coarseSize[i];

      Dune::FieldVector<ctype,dim> h;
      stream >>  input;
      for (int i=0; i<dim; i++)
        stream >> h[i];

      // the constructor takes the upper right corner...
      Dune::FieldVector<ctype,dim> length(h);
      for (int i=0; i<dim; i++)
        length[i] *= coarseSize[i];

      YLoadBalanceBackup<dim> lb(torus_dims);

      Grid* grid = new Dune::YaspGrid<dim>(length, coarseSize, periodic, overlap, comm, &lb);

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
    typedef typename Grid::Traits::CollectiveCommunication Comm;

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
      stream << "Torus structure: ";
      for (int i=0; i<dim; i++)
        stream << grid.torus().dims(i) << " ";
      stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: " << (grid.getRefineOption() ? "1" : "0") << std::endl;
      stream << "Coarse Size: ";
      for (int i=0; i<dim; i++)
        stream << grid.levelSize(0,i) << " ";
      stream << std::endl;

      grid.begin()->coords.print(stream);
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore (const std::string &filename, Comm comm = Comm())
    {
      std::ostringstream filename_str;
      filename_str << filename;
      filename_str << comm.rank();
      std::ifstream file(filename_str.str());
      if( file )
        return restore(file, comm);
      else
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename_str.str() << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore (std::istream &stream, Comm comm = Comm())
    {
      std::string input;

      Dune::array<int,dim> torus_dims;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> torus_dims[i];

      int refinement;
      stream >> input >> input;
      stream >> refinement;

      std::bitset<dim> periodic;
      bool b;
      stream >> input;
      for (int i=0; i<dim; i++)
      {
        stream >> b;
        periodic[i] = b;
      }

      int overlap;
      stream >> input;
      stream >> overlap;

      bool physicalOverlapSize;
      stream >> input;
      stream >> physicalOverlapSize;

      Dune::array<int,dim> coarseSize;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> coarseSize[i];

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

      YLoadBalanceBackup<dim> lb(torus_dims);
      Grid* grid = new Grid(coords, periodic, overlap, comm, coarseSize, &lb);

      grid->refineOptions(physicalOverlapSize);
      grid->globalRefine(refinement);

      return grid;
    }
  };
} // namespace Dune

#endif // #ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
