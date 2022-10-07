// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH

#include <cstdint>
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdint>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>

#include <dune/grid/io/file/vtk/streams.hh>
#include <dune/grid/io/file/vtk/common.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Data array writers for the VTKWriter

    This file contains classes to help writing data in the difeerent VTK
    output modes
 */

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! base class for data array writers
    /**
     * This is the base class; for an actual implementation look at
     * VTKAsciiDataArrayWriter, VTKBinaryDataArrayWriter, or
     * VTKBinaryAppendedDataArrayWriter.
     *
     * To create an actual DataArrayWriter, one would usually use an object of
     * class DataArrayWriterFactory.
     *
     * In the constructor, the actual writer implementation will prepare the
     * stream for writing the data.  This usually means writing a
     * "<DataArray>" header to the stream.  In the write() method it will
     * write a data element to the stream (this is not true for binaryappended
     * mode however: in this mode write just counts the number of bytes that
     * would be written, the actual writing has to happen later independent of
     * the writer).  Finally, in the destructor, the stream is put back in a
     * sane state.  That usually means writing something line "</DataArray>".
     */
    class DataArrayWriter
    {
    public:
      //! \brief construct a data array writer
      /**
       * \brief construct a data array writer
       * \param prec the precision type with which to write the data
       */
      DataArrayWriter(Precision _prec)
      : prec(_prec)
      {}

      //! write one element of data
      template<class T>
      void write(T data)
      {
        switch(prec)
        {
          case Precision::float32:
            writeFloat32(data); break;
          case Precision::float64:
            writeFloat64(data); break;
          case Precision::uint32:
            writeUInt32(data); break;
          case Precision::uint8:
            writeUInt8(data); break;
          case Precision::int32:
            writeInt32(data); break;
          default:
            DUNE_THROW(Dune::NotImplemented, "Unknown precision type");
        }
      }

      //! whether calls to write may be skipped
      virtual bool writeIsNoop() const { return false; }
      //! virtual destructor
      virtual ~DataArrayWriter () {}

    private:
      //! write one data element as float
      virtual void writeFloat32 (float data) = 0;
      //! write one data element as double
      virtual void writeFloat64 (double data) = 0;
      //! write one data element as integer
      virtual void writeInt32 (std::int32_t data) = 0;
      //! write one data element as unsigned integer
      virtual void writeUInt8 (std::uint8_t data) = 0;
      //! write one data element as unsigned integer
      virtual void writeUInt32 (std::uint32_t data) = 0;

      Precision prec;
    };

    //! a streaming writer for data array tags, uses ASCII inline format
    class AsciiDataArrayWriter : public DataArrayWriter
    {
    public:
      //! make a new data array writer
      /**
       * \param theStream Stream to write to.
       * \param name      Name of array to write.
       * \param ncomps    Number of components of the array.
       * \param indent_   Indentation to use.  This is use as-is for the
       *                  header and trailer lines, but increase by one level
       *                  for the actual data.
       */
      AsciiDataArrayWriter(std::ostream& theStream, std::string name,
                           int ncomps, const Indent& indent_, Precision prec_)
        : DataArrayWriter(prec_), s(theStream), counter(0), numPerLine(12), indent(indent_)
      {
        s << indent << "<DataArray type=\"" << toString(prec_) << "\" "
          << "Name=\"" << name << "\" ";
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"ascii\">\n";
        ++indent;
      }

      //! finish output; writes end tag
      ~AsciiDataArrayWriter ()
      {
        if (counter%numPerLine!=0) s << "\n";
        --indent;
        s << indent << "</DataArray>\n";
      }

    private:
      //! write one double data element to output stream
      void writeFloat64 (double data) final
      { write_float(data); }
      //! write one float data element to output stream
      void writeFloat32 (float data) final
      { write_float(data); }
      //! write one int data element to output stream
      void writeInt32 (std::int32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt32 (std::uint32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt8 (std::uint8_t data) final
      { write_(data); }

      template<class T>
      void write_(T data)
      {
        typedef typename PrintType<T>::Type PT;
        if(counter%numPerLine==0) s << indent;
        else s << " ";
        const auto original_precision = std::cout.precision();
        s << std::setprecision(std::numeric_limits<PT>::digits10) << (PT) data;
        std::cout.precision(original_precision);
        counter++;
        if (counter%numPerLine==0) s << "\n";
      }

      template<class T>
      void write_float(T data)
      {
        typedef typename PrintType<T>::Type PT;
        if(counter%numPerLine==0) s << indent;
        else s << " ";
        PT out_data = (PT) data;
        if (std::fpclassify(out_data) == FP_SUBNORMAL)
        {
          // truncate denormalized data to 0 to avoid Paraview segfaults on macOS
          out_data = 0;
        }
        const auto original_precision = std::cout.precision();
        s << std::setprecision(std::numeric_limits<PT>::digits10) << out_data;
        std::cout.precision(original_precision);
        counter++;
        if (counter%numPerLine==0) s << "\n";
      }

      std::ostream& s;
      int counter;
      int numPerLine;
      Indent indent;
    };

    //! a streaming writer for data array tags, uses binary inline format
    class BinaryDataArrayWriter : public DataArrayWriter
    {
    public:
      //! make a new data array writer
      /**
       * \param theStream Stream to write to.
       * \param name      Name of array to write.
       * \param ncomps    Number of components of the array.
       * \param nitems    Number of cells for cell data/Number of vertices for
       *                  point data.
       * \param indent_   Indentation to use.  This is use as-is for the
       *                  header and trailer lines, but increase by one level
       *                  for the actual data.
       */
      BinaryDataArrayWriter(std::ostream& theStream, std::string name,
                            int ncomps, int nitems, const Indent& indent_, Precision prec_)
        : DataArrayWriter(prec_), s(theStream), b64(theStream), indent(indent_)
      {
        s << indent << "<DataArray type=\"" << toString(prec_) << "\" "
          << "Name=\"" << name << "\" ";
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"binary\">\n";

        // write indentation for the data chunk
        s << indent+1;
        // store size, needs to be exactly 32 bit
        std::uint32_t size = ncomps*nitems*typeSize(prec_);
        b64.write(size);
        b64.flush();
      }

      //! finish output; writes end tag
      ~BinaryDataArrayWriter ()
      {
        b64.flush();
        // append newline to written data
        s << "\n";
        s << indent << "</DataArray>\n";
        s.flush();
      }

    private:
      //! write one double data element to output stream
      void writeFloat64 (double data) final
      { write_(data); }
      //! write one float data element to output stream
      void writeFloat32 (float data) final
      { write_(data); }
      //! write one int data element to output stream
      void writeInt32 (std::int32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt32 (std::uint32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt8 (std::uint8_t data) final
      { write_(data); }

      //! write one data element to output stream
      template<class T>
      void write_(T data)
      {
        b64.write(data);
      }

      std::ostream& s;
      Base64Stream b64;
      const Indent& indent;
    };

    //! a streaming writer for data array tags, uses appended raw format
    class AppendedRawDataArrayWriter : public DataArrayWriter
    {
    public:
      //! make a new data array writer
      /**
       * \param s         Stream to write to.
       * \param name      Name of array to write.
       * \param ncomps    Number of components of the array.
       * \param nitems    Number of cells for cell data/Number of vertices for
       *                  point data.
       * \param offset    Byte count variable: this is incremented by one for
       *                  each byte which has to written to the appended data
       *                  section later.
       * \param indent    Indentation to use.  This is uses as-is for the
       *                  header line.
       */
      AppendedRawDataArrayWriter(std::ostream& s, std::string name,
                                 int ncomps, unsigned nitems, unsigned& offset,
                                 const Indent& indent, Precision prec_)
      : DataArrayWriter(prec_)
      {
        s << indent << "<DataArray type=\"" << toString(prec_) << "\" "
          << "Name=\"" << name << "\" ";
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"appended\" offset=\""<< offset << "\" />\n";
        offset += 4; // header
        offset += ncomps*nitems*typeSize(prec_);
      }

      //! whether calls to write may be skipped
      bool writeIsNoop() const { return true; }

    private:
      //! write one data element to output stream (noop)
      void writeFloat64 (double) final {}
      void writeFloat32 (float) final {}
      void writeInt32 (std::int32_t) final {}
      void writeUInt32 (std::uint32_t) final {}
      void writeUInt8 (std::uint8_t) final {}
    };

    //! a streaming writer for data array tags, uses appended base64 format
    class AppendedBase64DataArrayWriter : public DataArrayWriter
    {
    public:
      //! make a new data array writer
      /**
       * \param s         Stream to write to.
       * \param name      Name of array to write.
       * \param ncomps    Number of components of the array.
       * \param nitems    Number of cells for cell data/Number of vertices for
       *                  point data.
       * \param offset    Byte count variable: this is incremented by one for
       *                  each base64 char which has to written to the
       *                  appended data section later.
       * \param indent    Indentation to use.  This is uses as-is for the
       *                  header line.
       */
      AppendedBase64DataArrayWriter(std::ostream& s, std::string name,
                                    int ncomps, unsigned nitems,
                                    unsigned& offset, const Indent& indent, Precision prec_)
      : DataArrayWriter(prec_)
      {
        s << indent << "<DataArray type=\"" << toString(prec_) << "\" "
          << "Name=\"" << name << "\" ";
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"appended\" offset=\""<< offset << "\" />\n";
        offset += 8; // header
        std::size_t bytes = ncomps*nitems*typeSize(prec_);
        offset += bytes/3*4;
        if(bytes%3 != 0)
          offset += 4;
      }

      //! whether calls to write may be skipped
      bool writeIsNoop() const { return true; }

    private:
      //! write one data element to output stream (noop)
      void writeFloat64 (double) final {}
      void writeFloat32 (float) final {}
      void writeInt32 (std::int32_t) final {}
      void writeUInt32 (std::uint32_t) final {}
      void writeUInt8 (std::uint8_t) final {}
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Naked ArrayWriters for the appended section
    //

    //! a streaming writer for appended data array tags, uses base64 format
    class NakedBase64DataArrayWriter : public DataArrayWriter
    {
    public:
      //! make a new data array writer
      /**
       * \param theStream Stream to write to.
       * \param ncomps    Number of components of the array.
       * \param nitems    Number of cells for cell data/Number of vertices for
       *                  point data.
       */
      NakedBase64DataArrayWriter(std::ostream& theStream, int ncomps,
                                 int nitems, Precision prec_)
        : DataArrayWriter(prec_), b64(theStream)
      {
        // store size
        std::uint32_t size = ncomps*nitems*typeSize(prec_);
        b64.write(size);
        b64.flush();
      }

    private:
      //! write one double data element to output stream
      void writeFloat64 (double data) final
      { write_(data); }
      //! write one float data element to output stream
      void writeFloat32 (float data) final
      { write_(data); }
      //! write one int data element to output stream
      void writeInt32 (std::int32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt32 (std::uint32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt8 (std::uint8_t data) final
      { write_(data); }

      //! write one data element to output stream
      template<class T>
      void write_(T data)
      {
          b64.write(data);
      }

      Base64Stream b64;
    };

    //! a streaming writer for appended data arrays, uses raw format
    class NakedRawDataArrayWriter : public DataArrayWriter
    {
      RawStream s;

    public:
      //! make a new data array writer
      /**
       * \param theStream Stream to write to.
       * \param ncomps    Number of components of the array.
       * \param nitems    Number of cells for cell data/Number of vertices for
       *                  point data.
       */
      NakedRawDataArrayWriter(std::ostream& theStream, int ncomps,
                              int nitems, Precision prec_)
        : DataArrayWriter(prec_), s(theStream)
      {
        s.write((unsigned int)(ncomps*nitems*typeSize(prec_)));
      }

    private:
      //! write one double data element to output stream
      void writeFloat64 (double data) final
      { write_(data); }
      //! write one float data element to output stream
      void writeFloat32 (float data) final
      { write_(data); }
      //! write one int data element to output stream
      void writeInt32 (std::int32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt32 (std::uint32_t data) final
      { write_(data); }
      //! write one unsigned int data element to output stream
      void writeUInt8 (std::uint8_t data) final
      { write_(data); }

      //! write one data element to output stream
      template<class T>
      void write_(T data)
      {
        s.write(data);
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Factory
    //

    //! a factory for DataArrayWriters
    /**
     * Some types of DataArrayWriters need to communicate data sauch as byte
     * counts between different instances which write to the same stream.
     * This factory will manage that data.
     */
    class DataArrayWriterFactory {
      enum Phase { main, appended };

      OutputType type;
      std::ostream& stream;
      unsigned offset;
      //! whether we are in the main or in the appended section writing phase
      Phase phase;

    public:
      //! create a DataArrayWriterFactory
      /**
       * \param type_   Type of DataArrayWriters to create
       * \param stream_ The stream that the DataArrayWriters will write to.
       *
       * Better avoid having multiple active factories on the same stream at
       * the same time.  Having an inactive factory (one whose make() method
       * is not called anymore before destruction) around at the same time as
       * an active one should be OK however.
       */
      inline DataArrayWriterFactory(OutputType type_, std::ostream& stream_)
        : type(type_), stream(stream_), offset(0), phase(main)
      { }

      //! signal start of the appended section
      /**
       * This method should be called after the main section has been written,
       * but before the appended section has been started.  After this method
       * has been called, a call to make will produce DataArrayWriters
       * suitable for the appended section.  The return value of this method
       * signals whether a appended section should be written at all: true
       * means write an appended section, false means don't write an appended
       * section.  If an appended section is not needed, the method make() may
       * not be called after a call to this method.
       */
      inline bool beginAppended() {
        phase = appended;
        switch(type) {
        case ascii :          return false;
        case base64 :         return false;
        case appendedraw :    return true;
        case appendedbase64 : return true;
        }
        DUNE_THROW(IOError, "Dune::VTK::DataArrayWriter: unsupported "
                   "OutputType " << type);
      }

      //! query encoding string for appended data
      const std::string& appendedEncoding() const {
        static const std::string rawString = "raw";
        static const std::string base64String = "base64";

        switch(type) {
        case ascii :
        case base64 :
          DUNE_THROW(IOError, "DataArrayWriterFactory::appendedEncoding(): No "
                     "appended encoding for OutputType " << type);
        case appendedraw :    return rawString;
        case appendedbase64 : return base64String;
        }
        DUNE_THROW(IOError, "DataArrayWriterFactory::appendedEncoding(): "
                   "unsupported OutputType " << type);
      }

      //! create a DataArrayWriter
      /**
       * \tparam T Type of the data to write.
       *
       * \param name   Name of the array to write.
       * \param ncomps Number of components of the vectors in to array.
       * \param nitems Number of vectors in the array.
       * \param indent Indentation to use.  This is use as-is for the header
       *               and trailer lines, but increase by one level for the
       *               actual data.
       * \param prec   the precision type of the output
       *
       * The should never be more than one DataArrayWriter on the same stream
       * around.  The returned object should be freed with delete.
       */
      DataArrayWriter* make(const std::string& name, unsigned ncomps,
                            unsigned nitems, const Indent& indent,
                            Precision prec)
      {
        switch(phase) {
        case main :
          switch(type) {
          case ascii :
            return new AsciiDataArrayWriter(stream, name, ncomps, indent, prec);
          case base64 :
            return new BinaryDataArrayWriter(stream, name, ncomps, nitems,
                                             indent, prec);
          case appendedraw :
            return new AppendedRawDataArrayWriter(stream, name, ncomps,
                                                  nitems, offset, indent, prec);
          case appendedbase64 :
            return new AppendedBase64DataArrayWriter(stream, name, ncomps,
                                                     nitems, offset,
                                                     indent, prec);
          }
          break;
        case appended :
          switch(type) {
          case ascii :
          case base64 :
            break; // invlid in appended mode
          case appendedraw :
            return new NakedRawDataArrayWriter(stream, ncomps, nitems, prec);
          case appendedbase64 :
            return new NakedBase64DataArrayWriter(stream, ncomps, nitems, prec);
          }
          break;
        }
        DUNE_THROW(IOError, "Dune::VTK::DataArrayWriter: unsupported "
                   "OutputType " << type << " in phase " << phase);
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
