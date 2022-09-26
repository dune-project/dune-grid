// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_PVTUWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_PVTUWRITER_HH

#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>

#include <dune/grid/io/file/vtk/common.hh>

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
         PVTUWriter writer(std::cout, polyData);

         // start the main section
         writer.beginMain();

         // dump cell data (optional)
         writer.beginCellData();
         for(each cell data field)
           writer.addArray(field.name, field.ncomps, precision);
         writer.endCellData();

         // dump point data (optional)
         writer.beginPointData();
         for(each point data field)
           writer.addArray(field.name, field.ncomps, precision);
         writer.endPointData();

         // dump point coordinates
         writer.beginPoints();
         writer.addArray("Coordinates", 3, precision);
         writer.endPoints();

         for(each serial piece)
           writer.addPiece(piece.filename);

         // finish main section
         writer.endMain();

         // end scope so the destructor gets called and the closing tag is written
       }
       \endcode
     */
    class PVTUWriter {
      std::ostream& stream;

      std::string fileType;

      Indent indent;

    public:
      //! create a PVTUWriter object
      /**
       * \param stream_    Stream to write to.
       * \param fileType_  Whether to write PolyData (1D) or UnstructuredGrid
       *                   (nD) format.
       *
       * Create object and write header.
       */
      inline PVTUWriter(std::ostream& stream_, FileType fileType_)
        : stream(stream_)
      {
        switch(fileType_) {
        case polyData :
          fileType = "PPolyData";
          break;
        case unstructuredGrid :
          fileType = "PUnstructuredGrid";
          break;
        default :
          DUNE_THROW(IOError, "PVTUWriter: Unknown fileType: " << fileType_);
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
      inline ~PVTUWriter() {
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
        stream << indent << "<PPointData";
        if(scalars != "") stream << " Scalars=\"" << scalars << "\"";
        if(vectors != "") stream << " Vectors=\"" << vectors << "\"";
        stream << ">\n";
        ++indent;
      }
      //! finish PointData section
      inline void endPointData() {
        --indent;
        stream << indent << "</PPointData>\n";
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
        stream << indent << "<PCellData";
        if(scalars != "") stream << " Scalars=\"" << scalars << "\"";
        if(vectors != "") stream << " Vectors=\"" << vectors << "\"";
        stream << ">\n";
        ++indent;
      }
      //! finish CellData section
      inline void endCellData() {
        --indent;
        stream << indent << "</PCellData>\n";
      }

      //! start section for the point coordinates
      /**
       * Between the call to this method an the following call to the
       * endPoints(), there must be a single field written.  The name must be
       * "Coordinates" and it must have 3 components.
       */
      inline void beginPoints() {
        stream << indent << "<PPoints>\n";
        ++indent;
      }
      //! finish section for the point coordinates
      inline void endPoints() {
        --indent;
        stream << indent << "</PPoints>\n";
      }

      //! start the main PPolyData/PUnstructuredGrid section
      /**
       * \param ghostLevel Set the GhostLevel attribute
       *
       * Between the call to this method and to endMain(), there should be
       * calls to add the actual data:
       * <ul>
       * <li> (optional) beginCellData()/endCellData(),
       * <li> (optional) beginPointData()/endPointData(),
       * <li> beginPoints()/endPoints(),
       * <li> one or more calls to addPiece()
       * </ul>
       */
      inline void beginMain(unsigned ghostLevel = 0) {
        stream << indent << "<" << fileType
               << " GhostLevel=\"" << ghostLevel << "\">\n";
        ++indent;
      }
      //! finish the main PolyData/UnstructuredGrid section
      inline void endMain() {
        --indent;
        stream << indent << "</" << fileType << ">\n";
      }

      //! Add an array to the output file
      /**
       * \tparam T The datatype of the array.
       * \param name   Name of the array.
       * \param ncomps Number of components in each vector of the array.
       * \param prec   The output precision of the data array
       */
      void addArray(const std::string& name, unsigned ncomps, Precision prec) {
        stream << indent << "<PDataArray"
               << " type=\"" << toString(prec) << "\""
               << " Name=\"" << name << "\""
               << " NumberOfComponents=\"" << ncomps << "\"/>\n";
      }

      //! Add a serial piece to the output file
      inline void addPiece(const std::string& filename) {
        stream << indent << "<Piece "
               << " Source=\"" << filename << "\"/>\n";
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_PVTUWRITER_HH
