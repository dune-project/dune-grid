// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_BASICWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_BASICWRITER_HH

#include <fstream>
#include <iomanip>
#include <iterator>
#include <list>
#include <memory>
#include <sstream>
#include <string>

#include <dune/common/parallel/mpiguard.hh>
#include <dune/common/path.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/functionwriter.hh>
#include <dune/grid/io/file/vtk/pvtuwriter.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    template<typename IteratorFactory>
    class BasicWriter {
      typedef typename IteratorFactory::CellIterator CellIterator;
      typedef typename IteratorFactory::CornerIterator CornerIterator;
      typedef typename IteratorFactory::PointIterator PointIterator;

      typedef typename IteratorFactory::Cell Cell;

    public:
      typedef FunctionWriterBase<Cell> FunctionWriter;

    private:
      typedef std::list<std::shared_ptr<FunctionWriter> > WriterList;
      typedef typename WriterList::const_iterator WIterator;

      typedef typename Cell::Geometry::ctype ctype;
      static const unsigned celldim = Cell::mydimension;
      typedef ReferenceElements<ctype, celldim> Refelems;

      static const FileType fileType = celldim == 1
                                       ? polyData : unstructuredGrid;

      const IteratorFactory& factory;

      WriterList cellData;
      WriterList pointData;

      CoordinatesWriter<Cell> coords;
      typename IteratorFactory::ConnectivityWriter connectivity;
      OffsetsWriter<Cell> offsets;
      TypesWriter<Cell> types;

    public:
      BasicWriter(const IteratorFactory& factory_)
        : factory(factory_), connectivity(factory.makeConnectivity())
      { }

      //////////////////////////////////////////////////////////////////////
      //
      //  Methods for adding data
      //

      void addCellData(const std::shared_ptr<FunctionWriter>& writer) {
        cellData.push_back(writer);
      }

      void addPointData(const std::shared_ptr<FunctionWriter>& writer) {
        pointData.push_back(writer);
      }

      void clear() {
        cellData.clear();
        pointData.clear();
      }

    protected:
      //////////////////////////////////////////////////////////////////////
      //
      //  Methods for writing single functions
      //

      void writeCellFunction(VTUWriter& vtuWriter,
                             FunctionWriter& functionWriter,
                             unsigned ncells) const
      {
        if(functionWriter.beginWrite(vtuWriter, ncells)) {
          const CellIterator& cellend = factory.endCells();
          for(CellIterator cellit = factory.beginCells(); cellit != cellend;
              ++cellit)
            functionWriter.write(*cellit, Refelems::general(cellit->type()).
                                 position(0,0));
        }
        functionWriter.endWrite();
      }

      void writePointFunction(VTUWriter& vtuWriter,
                              FunctionWriter& functionWriter,
                              unsigned npoints) const
      {
        if(functionWriter.beginWrite(vtuWriter, npoints)) {
          const PointIterator& pend = factory.endPoints();
          for(PointIterator pit = factory.beginPoints(); pit != pend; ++pit)
            functionWriter.write(pit->cell(), pit->duneIndex());
        }
        functionWriter.endWrite();
      }

      void writeCornerFunction(VTUWriter& vtuWriter,
                               FunctionWriter& functionWriter,
                               unsigned ncorners) const
      {
        if(functionWriter.beginWrite(vtuWriter, ncorners)) {
          const CornerIterator& cend = factory.endCorners();
          for(CornerIterator cit = factory.beginCorners(); cit != cend; ++cit)
            functionWriter.write(cit->cell(), cit->duneIndex());
        }
        functionWriter.endWrite();
      }

      //////////////////////////////////////////////////////////////////////
      //
      //  Methods for writing whole sections
      //

      static std::string getFirstScalar(const WriterList& data) {
        const WIterator& wend = data.end();
        for(WIterator wit = data.begin(); wit != wend; ++wit)
          if((*wit)->ncomps() == 1)
            return (*wit)->name();
        return "";
      }

      static std::string getFirstVector(const WriterList& data) {
        const WIterator& wend = data.end();
        for(WIterator wit = data.begin(); wit != wend; ++wit)
          if((*wit)->ncomps() == 3)
            return (*wit)->name();
        return "";
      }

      void writeCellData(VTUWriter& vtuWriter, unsigned ncells) const {
        if(cellData.empty()) return;

        vtuWriter.beginCellData(getFirstScalar(cellData),
                                getFirstVector(cellData));
        const WIterator& wend = cellData.end();
        for(WIterator wit = cellData.begin(); wit != wend; ++wit)
          writeCellFunction(vtuWriter, **wit, ncells);
        vtuWriter.endCellData();
      }

      void writePointData(VTUWriter& vtuWriter, unsigned npoints) const {
        if(pointData.empty()) return;

        vtuWriter.beginPointData(getFirstScalar(pointData),
                                 getFirstVector(pointData));
        const WIterator& wend = pointData.end();
        for(WIterator wit = pointData.begin(); wit != wend; ++wit)
          writePointFunction(vtuWriter, **wit, npoints);
        vtuWriter.endPointData();
      }

      void writeGrid(VTUWriter& vtuWriter, unsigned ncells, unsigned npoints,
                     unsigned ncorners) {
        vtuWriter.beginPoints();
        writePointFunction(vtuWriter, coords, npoints);
        vtuWriter.endPoints();

        vtuWriter.beginCells();
        writeCornerFunction(vtuWriter, connectivity, ncorners);
        writeCellFunction(vtuWriter, offsets, ncells);
        if(fileType != polyData)
          writeCellFunction(vtuWriter, types, ncells);
        vtuWriter.endCells();
      }

      void writeAll(VTUWriter& vtuWriter, unsigned ncells, unsigned npoints,
                    unsigned ncorners) {
        writeCellData(vtuWriter, ncells);
        writePointData(vtuWriter, npoints);
        writeGrid(vtuWriter, ncells, npoints, ncorners);
      }

    public:
      void writePiece(const std::string& filename, OutputType outputType) {
        std::ofstream stream;
        stream.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                          std::ios_base::eofbit);
        stream.open(filename.c_str(), std::ios::binary);

        VTUWriter vtuWriter(stream, outputType, fileType);

        unsigned ncells = std::distance(factory.beginCells(),
                                        factory.endCells());
        unsigned npoints = std::distance(factory.beginPoints(),
                                         factory.endPoints());
        unsigned ncorners = std::distance(factory.beginCorners(),
                                          factory.endCorners());

        vtuWriter.beginMain(ncells, npoints);
        writeAll(vtuWriter, ncells, npoints, ncorners);
        vtuWriter.endMain();

        if(vtuWriter.beginAppended())
          writeAll(vtuWriter, ncells, npoints, ncorners);
        vtuWriter.endAppended();

      }

      //! write header file in parallel case to stream
      /**
       * Writes a .pvtu/.pvtp file for a collection of concurrently written
       * .vtu/.vtp files.
       *
       * \param name      Name of file to write contents to,

       * \param piecename Base name of the pieces.  Should not contain a
       *                  directory part or filename extension.
       * \param piecepath Directory part of the pieces.  Since paraview does
       *                  not support absolute paths in parallel collection
       *                  files, this should be a path relative to the
       *                  directory the collection file resides in.  A
       *                  trailing '/' is optional, and an empty value "" is
       *                  equivalent to the value "." except it will look
       *                  nicer in the collection file.
       */
      void writeCollection(const std::string name,
                           const std::string& piecename,
                           const std::string& piecepath)
      {
        std::ofstream stream;
        stream.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                          std::ios_base::eofbit);
        stream.open(name.c_str(), std::ios::binary);

        PVTUWriter writer(stream, fileType);

        writer.beginMain();

        // PPointData
        writer.beginPointData(getFirstScalar(pointData),
                              getFirstVector(pointData));
        for(WIterator it=pointData.begin(); it!=pointData.end(); ++it)
          (*it)->addArray(writer);
        writer.endPointData();

        // PCellData
        writer.beginCellData(getFirstScalar(cellData),
                             getFirstVector(cellData));
        for(WIterator it=cellData.begin(); it!=cellData.end(); ++it)
          (*it)->addArray(writer);
        writer.endCellData();

        // PPoints
        writer.beginPoints();
        coords.addArray(writer);
        writer.endPoints();

        // Pieces
        for( int i = 0; i < factory.comm().size(); ++i )
          writer.addPiece(getParallelPieceName(piecename, piecepath, i));

        writer.endMain();
      }

      //////////////////////////////////////////////////////////////////////
      //
      //  Filename generators
      //

      //! return name of a parallel piece file
      /**
       * \param name Base name of the VTK output.  This should be without any
       *             directory parts and without a filename extension.
       * \param path Directory part of the resulting piece name.  May be
       *             empty, in which case the resulting name will not have a
       *             directory part.  If non-empty, may or may not have a
       *             trailing '/'.  If a trailing slash is missing, one is
       *             appended implicitly.
       * \param rank Rank of the process to generate a piece name for.
       */
      std::string getParallelPieceName(const std::string& name,
                                       const std::string& path, int rank) const
      {
        std::ostringstream s;
        if(path.size() > 0) {
          s << path;
          if(path[path.size()-1] != '/')
            s << '/';
        }
        s << 's' << std::setw(4) << std::setfill('0') << factory.comm().size()
          << ':';
        s << 'p' << std::setw(4) << std::setfill('0') << rank << ':';
        s << name;
        switch(fileType) {
        case polyData :         s << ".vtp"; break;
        case unstructuredGrid : s << ".vtu"; break;
        }
        return s.str();
      }

      //! return name of a parallel header file
      /**
       * \param name Base name of the VTK output.  This should be without any
       *             directory parts and without a filename extension.
       * \param path Directory part of the resulting header name.  May be
       *             empty, in which case the resulting name will not have a
       *             directory part.  If non-empty, may or may not have a
       *             trailing '/'.  If a trailing slash is missing, one is
       *             appended implicitly.
       */
      std::string getParallelHeaderName(const std::string& name,
                                        const std::string& path) const
      {
        std::ostringstream s;
        if(path.size() > 0) {
          s << path;
          if(path[path.size()-1] != '/')
            s << '/';
        }
        s << 's' << std::setw(4) << std::setfill('0') << factory.comm().size()
          << ':';
        s << name;
        switch(fileType) {
        case polyData :         s << ".pvtp"; break;
        case unstructuredGrid : s << ".pvtu"; break;
        }
        return s.str();
      }

      //! return name of a serial piece file
      /**
       * This is similar to getParallelPieceName, but skips the prefixes for
       * commSize ("s####:") and commRank ("p####:").
       *
       * \param name Base name of the VTK output.  This should be without any
       *             directory parts and without a filename extension.
       * \param path Directory part of the resulting piece name.  May be
       *             empty, in which case the resulting name will not have a
       *             directory part.  If non-empty, may or may not have a
       *             trailing '/'.  If a trailing slash is missing, one is
       *             appended implicitly.
       */
      std::string getSerialPieceName(const std::string& name,
                                     const std::string& path) const
      {
        switch(fileType) {
        case polyData :         return concatPaths(path, name+".vtp");
        case unstructuredGrid : return concatPaths(path, name+".vtu");
        }
        return concatPaths(path, name); // unknown fileType
      }

      //////////////////////////////////////////////////////////////////////
      //
      //  User interface functions for writing
      //

      //! write output; interface might change later
      /**
       * \param name       Base name of the output files.  This should not
       *                   contain any directory part and not filename
       *                   extensions.  It will be used both for each processes
       *                   piece as well as the parallel collection file.
       * \param path       Directory where to put the parallel collection
       *                   (.pvtu/.pvtp) file.  If it is relative, it is taken
       *                   realtive to the current directory.
       * \param extendpath Directory where to put the piece file (.vtu/.vtp) of
       *                   this process.  If it is relative, it is taken
       *                   relative to the directory denoted by path.
       * \param outputType How to encode the data in the file.
       *
       * \note Currently, extendpath may not be absolute unless path is
       *       absolute, because that would require the value of the current
       *       directory.
       *
       * \throw NotImplemented Extendpath is absolute but path is relative.
       * \throw IOError        Failed to open a file.
       * \throw MPIGuardError  An exception was thrown during this method in
       *                       one of the other processes.
       */
      std::string pwrite(const std::string& name, const std::string& path,
                         const std::string& extendpath, OutputType outputType)
      {
        MPIGuard guard(factory.comm());

        // do some magic because paraview can only cope with relative paths to
        // piece files
        std::ofstream file;
        file.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                        std::ios_base::eofbit);
        std::string piecepath = concatPaths(path, extendpath);
        std::string relpiecepath = relativePath(path, piecepath);

        // write this processes .vtu/.vtp piece file
        std::string fullname = getParallelPieceName(name, piecepath,
                                                    factory.comm().rank());
        writePiece(fullname, outputType);

        // if we are rank 0, write .pvtu/.pvtp parallel header
        fullname = getParallelHeaderName(name, path);
        if(factory.comm().rank() == 0)
          writeCollection(fullname, name, relpiecepath);

        guard.finalize();

        return fullname;
      }

      /** \brief write output (interface might change later)
       *
       *  This method can  be used in parallel as well  as in serial programs.
       *  For  serial runs  (commSize=1) it  chooses other  names  without the
       *  "s####:p####:" prefix  for the .vtu/.vtp files and  omits writing of
       *  the .pvtu/pvtp file however.  For parallel runs (commSize > 1) it is
       *  the same as a call to pwrite() with path="" and extendpath="".
       *
       *  \param name     Base name of the output files.  This should not
       *                  contain any directory part and no filename
       *                  extensions.
       *  \param outputType How to encode the data in the file.
       */
      std::string write(const std::string &name, OutputType outputType)
      {
        // in the parallel case, just use pwrite, it has all the necessary
        // stuff, so we don't need to reimplement it here.
        if(factory.comm().size() > 1)
          return pwrite(name, "", "", outputType);

        // generate filename for process data
        std::string pieceName = getSerialPieceName(name, "");

        writePiece(pieceName, outputType);

        return pieceName;
      }

    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_BASICWRITER_HH
