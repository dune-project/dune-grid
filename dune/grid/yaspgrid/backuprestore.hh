// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
#define DUNE_GRID_YASPGRID_BACKUPRESTORE_HH

//- system headers
#include <fstream>
#include <iterator>

//- Dune headers
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/yaspgrid.hh>

// bump this version number up if you introduce any changes
// to the outout format of the YaspGrid BackupRestoreFacility.
#define YASPGRID_BACKUPRESTORE_FORMAT_VERSION 2

namespace Dune
{

  template<class Coordinates>
  struct MaybeHaveOrigin
  {

    template<class S>
    static void writeOrigin(S& /* s */, const Coordinates& /* coord */)
    {}

    template<class S>
    static void readOrigin(S& /* s */, Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& /* coord */)
    {}

    template<typename... A>
    static typename Dune::YaspGrid<Coordinates::dimension,Coordinates>* createGrid(
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& /* lowerleft */, A... args)
    {
      return new Dune::YaspGrid<Coordinates::dimension,Coordinates>(args...);
    }
  };

  template<class ctype, int dim>
  struct MaybeHaveOrigin<Dune::EquidistantOffsetCoordinates<ctype, dim> >
  {
    typedef typename Dune::EquidistantOffsetCoordinates<ctype, dim> Coordinates;

    template<class S>
    static void writeOrigin(S& s, const Coordinates& coord)
    {
      s << "Origin: ";
      for (int i=0; i<dim; i++)
        s << coord.origin(i) << " ";
      s << std::endl;
    }

    template<class S>
    static void readOrigin(S& s, Dune::FieldVector<ctype, dim>& coord)
    {
      std::string token;
      s >> token;
      for (int i=0; i<dim; i++)
        s >> coord[i];
    }

    template<typename... A>
    static typename Dune::YaspGrid<Coordinates::dimension,Coordinates>* createGrid(
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& lowerleft,
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& extension, A... args)
    {
      Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension> upperright(lowerleft);
      upperright += extension;
      return new Dune::YaspGrid<Coordinates::dimension,Coordinates>(lowerleft, upperright, args...);
    }
  };

  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class Coordinates>
  struct BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >
  {
    // type of grid
    typedef typename Dune::YaspGrid<dim, Coordinates> Grid;
    typedef typename Grid::ctype ctype;
    typedef typename Grid::Traits::Communication Comm;

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
      stream << "YaspGrid BackupRestore Format Version: " << YASPGRID_BACKUPRESTORE_FORMAT_VERSION << std::endl;
      stream << "Torus structure: ";
      for (int i=0; i<dim; i++)
        stream << grid.torus().dims(i) << " ";
      stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: ";
      for (typename Grid::YGridLevelIterator i=std::next(grid.begin()); i != grid.end(); ++i)
        stream << (i->keepOverlap ? "1" : "0") << " ";
      stream << std::endl;
      stream << "Coarse Size: ";
      for (int i=0; i<dim; i++)
        stream << grid.levelSize(0,i) << " ";
      stream << std::endl;
      stream << "Meshsize: " ;
      for (int i=0; i<dim; i++)
        stream << grid.begin()->coords.meshsize(i,0) << " ";
      stream << std::endl;
      MaybeHaveOrigin<Coordinates>::writeOrigin(stream, grid.begin()->coords);
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

      int version;
      stream >> input >> input >> input >> input;
      stream >> version;
      if (version != YASPGRID_BACKUPRESTORE_FORMAT_VERSION)
        DUNE_THROW(Dune::Exception, "Your YaspGrid backup file is written in an outdated format!");

      std::array<int,dim> torus_dims;
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

      std::vector<bool> physicalOverlapSize;
      physicalOverlapSize.resize(refinement);
      stream >> input;
      for (int i=0; i<refinement; ++i)
      {
        stream >> b;
        physicalOverlapSize[i] = b;
      }

      std::array<int,dim> coarseSize;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> coarseSize[i];

      Dune::FieldVector<ctype,dim> h;
      stream >>  input;
      for (int i=0; i<dim; i++)
        stream >> h[i];

      Dune::FieldVector<ctype,dim> origin;
      MaybeHaveOrigin<Coordinates>::readOrigin(stream, origin);

      // the constructor takes the upper right corner...
      Dune::FieldVector<ctype,dim> length(h);
      for (int i=0; i<dim; i++)
        length[i] *= coarseSize[i];

      YaspFixedSizePartitioner<dim> lb(torus_dims);

      Grid* grid = MaybeHaveOrigin<Coordinates>::createGrid(origin, length, coarseSize, periodic, overlap, comm, &lb);

      for (int i=0; i<refinement; ++i)
      {
        grid->refineOptions(physicalOverlapSize[i]);
        grid->globalRefine(1);
      }

      return grid;
    }
  };

  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class ctype>
  struct BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,TensorProductCoordinates<ctype,dim> > Grid;
    typedef typename Grid::Traits::Communication Comm;

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
      stream << "YaspGrid BackupRestore Format Version: " << YASPGRID_BACKUPRESTORE_FORMAT_VERSION << std::endl;
      stream << "Torus structure: ";
      for (int i=0; i<dim; i++)
        stream << grid.torus().dims(i) << " ";
      stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
      stream << "Periodicity: ";
      for (int i=0; i<dim; i++)
        stream << (grid.isPeriodic(i) ? "1 " : "0 ");
      stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
      stream << "KeepPhysicalOverlap: ";
      for (typename Grid::YGridLevelIterator i=std::next(grid.begin()); i != grid.end(); ++i)
        stream << (i->keepOverlap ? "1" : "0") << " ";
      stream << std::endl;
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

      int version;
      stream >> input >> input >> input >> input;
      stream >> version;
      if (version != YASPGRID_BACKUPRESTORE_FORMAT_VERSION)
        DUNE_THROW(Dune::Exception, "Your YaspGrid backup file is written in an outdated format!");

      std::array<int,dim> torus_dims;
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

      std::vector<bool> physicalOverlapSize;
      physicalOverlapSize.resize(refinement);
      stream >> input;
      for (int i=0; i<refinement; ++i)
      {
        stream >> b;
        physicalOverlapSize[i] = b;
      }


      std::array<int,dim> coarseSize;
      stream >> input >> input;
      for (int i=0; i<dim; i++)
        stream >> coarseSize[i];

      std::array<std::vector<ctype>,dim> coords;
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

      YaspFixedSizePartitioner<dim> lb(torus_dims);
      Grid* grid = new Grid(coords, periodic, overlap, comm, coarseSize, &lb);

      for (int i=0; i<refinement; ++i)
      {
        grid->refineOptions(physicalOverlapSize[i]);
        grid->globalRefine(1);
      }

      return grid;
    }
  };
} // namespace Dune

#endif // #ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
