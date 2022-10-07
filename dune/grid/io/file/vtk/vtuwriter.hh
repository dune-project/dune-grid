// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_VTUWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_VTUWRITER_HH

#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>

namespace Dune {

  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! Dump a .vtu/.vtp files contents to a stream
    /**
     * This will help generating a .vtu/.vtp file.  Typical use is like this:
       \code
       {
         // create writer, writes begin tag
         VTUWriter writer(std::cout, appendedraw, VTUWriter::polyData);

         // write the main header
         writer.beginMain(ncells, nvertices);
         dumpEverything(writer);
         writer.endMain();

         // write the appended section, if required by the outputtype
         if(writer.beginAppended())
           dumpEverything(writer);
         writer.endAppended();

         // end scope so the destructor gets called and the closing tag is written
       }
       \endcode
     * The method dumpEverything() then looks something like this:
     * \code
       void dumpEverything(VTUWriter& writer)
       {
         // dump cell data (optional)
         writer.beginCellData();
         for(each cell data field) {
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter(field.name, field.ncomps, ncells, precision));
           // iterate over the points and write data for each
         }
         writer.endCellData();

         // dump point data (optional)
         writer.beginPointData();
         for(each point data field) {
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter(field.name, field.ncomps, npoints, precision));
           // iterate over the points and write data for each
         }
         writer.endPointData();

         // dump point coordinates
         writer.beginPoints();
         {
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter("Coordinates", 3, npoints, precision));
           // iterate over the points and write data for each
         }
         writer.endPoints();

         // dump cells
         writer.beginCells();
         { // connectivity
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter("connectivity", 1, ncorners, precision));
           // iterate over the cells and write data for each
         }
         { // connectivity
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter("offsets", 1, ncells, precision));
           // iterate over the cells and write data for each
         }
         if(fileType == unstructuredGrid) { // types
           shared_ptr<DataArrayWriter> arraywriter
             (writer.makeArrayWriter("types", 1, ncells, precision));
           // iterate over the cells and write data for each
         }
         writer.endCells();
       }
       \endcode
     */
    class VTUWriter {
    public:
      std::ostream& stream;
      enum Phase { main, appended } phase;

    private:
      DataArrayWriterFactory factory;
      Indent indent;

      std::string fileType;
      std::string cellName;

      bool doAppended;

    public:
      //! create a VTUWriter object
      /**
       * \param stream_    Stream to write to.
       * \param outputType How to encode data.
       * \param fileType_  Whether to write PolyData (1D) or UnstructuredGrid
       *                   (nD) format.
       *
       * Create object and write header.
       */
      inline VTUWriter(std::ostream& stream_, OutputType outputType,
                       FileType fileType_)
        : stream(stream_), factory(outputType, stream)
      {
        switch(fileType_) {
        case polyData :
          fileType = "PolyData";
          cellName = "Lines";
          break;
        case unstructuredGrid :
          fileType = "UnstructuredGrid";
          cellName = "Cells";
          break;
        default :
          DUNE_THROW(IOError, "VTUWriter: Unknown fileType: " << fileType_);
        }
        const std::string& byteOrder = getEndiannessString();

        stream << indent << "<?xml version=\"1.0\"?>\n";
        stream << indent << "<VTKFile"
               << " type=\"" << fileType << "\""
               << " version=\"0.1\""
               << " byte_order=\"" << byteOrder << "\">\n";
        ++indent;
      }

      //! write footer
      inline ~VTUWriter() {
        --indent;
        stream << indent << "</VTKFile>\n"
               << std::flush;
      }

      //! start PointData section
      /**
       * \param scalars Name of field to which should be marked as default
       *                scalars field.  If this is the empty string, don't set
       *                any default.
       * \param vectors Name of field to which should be marked as default
       *                vectors field.  If this is the empty string, don't set
       *                any default.
       *
       * If there are no PointData fields, the call to this function may be
       * skipped, together with the corresponding call to endPointData().
       */
      inline void beginPointData(const std::string& scalars = "",
                                 const std::string& vectors = "") {
        switch(phase) {
        case main :
          stream << indent << "<PointData";
          if(scalars != "") stream << " Scalars=\"" << scalars << "\"";
          if(vectors != "") stream << " Vectors=\"" << vectors << "\"";
          stream << ">\n";
          ++indent;
          break;
        case appended :
          break;
        }
      }
      //! finish PointData section
      inline void endPointData() {
        switch(phase) {
        case main :
          --indent;
          stream << indent << "</PointData>\n";
          break;
        case appended :
          break;
        }
      }

      //! start CellData section
      /**
       * \param scalars Name of field to which should be marked as default
       *                scalars field.  If this is the empty string, don't set
       *                any default.
       * \param vectors Name of field to which should be marked as default
       *                vectors field.  If this is the empty string, don't set
       *                any default.
       *
       * If there are no CellData fields, the call to this function may be
       * skipped, together with the corresponding call to endCellData().
       */
      inline void beginCellData(const std::string& scalars = "",
                                const std::string& vectors = "") {
        switch(phase) {
        case main :
          stream << indent << "<CellData";
          if(scalars != "") stream << " Scalars=\"" << scalars << "\"";
          if(vectors != "") stream << " Vectors=\"" << vectors << "\"";
          stream << ">\n";
          ++indent;
          break;
        case appended :
          break;
        }
      }
      //! finish CellData section
      inline void endCellData() {
        switch(phase) {
        case main :
          --indent;
          stream << indent << "</CellData>\n";
          break;
        case appended :
          break;
        }
      }

      //! start section for the point coordinates
      /**
       * Between the call to this method an the following call to the
       * endPoints(), there must be a single field written.  The name must be
       * "Coordinates", it must have 3 components, and the number of items
       * must be the number of points.
       */
      inline void beginPoints() {
        switch(phase) {
        case main :
          stream << indent << "<Points>\n";
          ++indent;
          break;
        case appended :
          break;
        }
      }
      //! finish section for the point coordinates
      inline void endPoints() {
        switch(phase) {
        case main :
          --indent;
          stream << indent << "</Points>\n";
          break;
        case appended :
          break;
        }
      }

      //! start section for the grid cells/PolyData lines
      /**
       * Between the call to this method an the following call to the
       * endCells(), there must be two or three fields written:
       * <ul>
       * <li>"connectivity" of type Int32 with 3 components, number of items
       *     is the number of corners (that may be different from number of
       *     vertices!)
       * <li>"offsets" of type Int32 with one component, number of items is
       *     number of cells.
       * <li>for UnstructuredGrid, "types" of type UInt8 with one component,
       *     number of items is number of cells.
       * </ul>
       */
      inline void beginCells() {
        switch(phase) {
        case main :
          stream << indent << "<" << cellName << ">\n";
          ++indent;
          break;
        case appended :
          break;
        }
      }
      //! start section for the grid cells/PolyData lines
      inline void endCells() {
        switch(phase) {
        case main :
          --indent;
          stream << indent << "</" << cellName << ">\n";
          break;
        case appended :
          break;
        }
      }

      //! start the main PolyData/UnstructuredGrid section
      /**
       * \param ncells  Number of cells/lines.
       * \param npoints Number of points.
       *
       * Between the call to this method and to endMain(), there should be
       * calls to dump the actual data:
       * <ul>
       * <li> (optional) beginCellData()/endCellData(),
       * <li> (optional) beginPointData()/endPointData(),
       * <li> beginPoints()/endPoints(),
       * <li> beginCells()/endCells(),
       * </ul>
       */
      inline void beginMain(unsigned ncells, unsigned npoints) {
        stream << indent << "<" << fileType << ">\n";
        ++indent;
        stream << indent << "<Piece"
               << " NumberOf" << cellName << "=\"" << ncells << "\""
               << " NumberOfPoints=\"" << npoints << "\">\n";
        ++indent;
        phase = main;
      }
      //! finish the main PolyData/UnstructuredGrid section
      inline void endMain() {
        --indent;
        stream << indent << "</Piece>\n";
        --indent;
        stream << indent << "</" << fileType << ">\n";
      }

      //! start the appended data section
      /**
       * \returns A value indicating whether the is an actual appended section
       *          required.
       *
       * If this function returns true, an appended section is actually
       * required.  In this case, between the call to this method and to
       * endAppended(), there should be literally the same calls (including
       * the same arguments) as between the calls to beginMain() and
       * endMain().  The only exception is, that if a DataArrayWriter in the
       * main section indicated that the calls to write could be skipped, this
       * is not necessarily true in the appended section also (you will have
       * to ask the DataArrayWriter again).
       *
       * If this function returns false, no appended section is required and a
       * call to endAppeded() should immediately follow the call to this
       * function.
       */
      inline bool beginAppended() {
        doAppended = factory.beginAppended();
        if(doAppended) {
          const std::string& encoding = factory.appendedEncoding();
          stream << indent << "<AppendedData"
                 << " encoding=\"" << encoding << "\">\n";
          ++indent;
          // mark begin of data
          stream << indent << "_";
        }
        phase = appended;
        return doAppended;
      }
      //! finish the appended data section
      inline void endAppended() {
        if(doAppended) {
          stream << "\n";
          --indent;
          stream << indent << "</AppendedData>\n";
        }
      }

      //! acquire a DataArrayWriter
      /**
       * \tparam T Type of the data to write.
       *
       * \param name   Name of the array to write.
       * \param ncomps Number of components of the vectors in the array.
       * \param nitems Number of vectors in the array (number of cells/number
       *               of points/number of corners).
       *
       * There should never be more than one DataArrayWriter created by the
       * same VTUWriter around.  The returned object should be freed with
       * delete.
       */
      DataArrayWriter* makeArrayWriter(const std::string& name,
                                       unsigned ncomps, unsigned nitems,
                                       Precision prec) {
        return factory.make(name, ncomps, nitems, indent, prec);
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_VTUWRITER_HH
