// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH

#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>

#include <dune/grid/io/file/vtk/streams.hh>
#include <dune/grid/io/file/vtk/common.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Data array wrriters for the VTKWriter

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
     * \tparam T Type of the data elements to write
     *
     * This is an abstract base class; for an actual implementation look at
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
    template<class T>
    class DataArrayWriter
    {
    public:
      //! write one data element
      virtual void write (T data) = 0;
      //! virtual destructor
      virtual ~DataArrayWriter () {}
    };

    //! a streaming writer for data array tags, uses ASCII inline format
    template<class T>
    class AsciiDataArrayWriter : public DataArrayWriter<T>
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
                           int ncomps, const Indent& indent_)
        : s(theStream), counter(0), numPerLine(12), indent(indent_)
      {
        TypeName<T> tn;
        s << indent << "<DataArray type=\"" << tn() << "\" "
          << "Name=\"" << name << "\" ";
        // vtk file format: a vector data always should have 3 comps (with 3rd
        // comp = 0 in 2D case)
        if (ncomps>3)
          DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                     "components");
        s << "NumberOfComponents=\"" << (ncomps>1 ? 3 : 1) << "\" ";
        s << "format=\"ascii\">\n";
        ++indent;
      }

      //! write one data element to output stream
      void write (T data)
      {
        typedef typename PrintType<T>::Type PT;
        if(counter%numPerLine==0) s << indent;
        else s << " ";
        s << (PT) data;
        counter++;
        if (counter%numPerLine==0) s << "\n";
      }

      //! finish output; writes end tag
      ~AsciiDataArrayWriter ()
      {
        if (counter%numPerLine!=0) s << "\n";
        --indent;
        s << indent << "</DataArray>\n";
      }

    private:
      std::ostream& s;
      int counter;
      int numPerLine;
      Indent indent;
    };

    //! a streaming writer for data array tags, uses binary inline format
    template<class T>
    class BinaryDataArrayWriter : public DataArrayWriter<T>
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
                            int ncomps, int nitems, const Indent& indent_)
        : s(theStream), b64(theStream), indent(indent_)
      {
        TypeName<T> tn;
        ncomps = (ncomps>1 ? 3 : 1);
        s << indent << "<DataArray type=\"" << tn() << "\" "
          << "Name=\"" << name << "\" ";
        // vtk file format: a vector data always should have 3 comps (with 3rd
        // comp = 0 in 2D case)
        if (ncomps>3)
          DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                     "components");
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"binary\">\n";

        // write indentation for the data chunk
        s << indent+1;
        // store size
        unsigned long int size = ncomps*nitems*sizeof(T);
        b64.write(size);
        b64.flush();
      }

      //! write one data element to output stream
      void write (T data)
      {
        b64.write(data);
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
      std::ostream& s;
      Base64Stream b64;
      const Indent& indent;
    };

    //! a streaming writer for data array tags, uses binary appended format
    template<class T>
    class BinaryAppendedDataArrayWriter : public DataArrayWriter<T>
    {
    public:
      //! make a new data array writer
      /**
       * \param theStream Stream to write to.
       * \param name      Name of array to write.
       * \param ncomps    Number of components of the array.
       * \param bc        Byte count variable: this is incremented by one for
       *                  each byte which has to beritte to the appended data
       *                  section later.
       * \param indent_   Indentation to use.  This is uses as-is for the
       *                  header line.
       */
      BinaryAppendedDataArrayWriter(std::ostream& theStream, std::string name,
                                    int ncomps, unsigned int& bc,
                                    const Indent& indent)
        : s(theStream),bytecount(bc)
      {
        TypeName<T> tn;
        s << indent << "<DataArray type=\"" << tn() << "\" "
          << "Name=\"" << name << "\" ";
        // vtk file format: a vector data always should have 3 comps(with 3rd
        // comp = 0 in 2D case)
        if (ncomps>3)
          DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                     "components");
        s << "NumberOfComponents=\"" << (ncomps>1 ? 3 : 1) << "\" ";
        s << "format=\"appended\" offset=\""<< bytecount << "\" />\n";
        bytecount += 4; // header
      }

      //! write one data element to output stream
      void write (T data)
      {
        bytecount += sizeof(T);
      }

    private:
      std::ostream& s;
      unsigned int& bytecount;
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Naked ArrayWriters for the appended section
    //

    //! a streaming writer for appended data array tags, uses base64 format
    template<class T>
    class NakedBase64DataArrayWriter : public DataArrayWriter<T>
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
                                 int nitems)
        : b64(theStream)
      {
        // store size
        unsigned long int size = ncomps*nitems*sizeof(T);
        b64.write(size);
        b64.flush();
      }

      //! write one data element to output stream
      void write (T data)
      {
        b64.write(data);
      }

    private:
      Base64Stream b64;
    };

    //! a streaming writer for appended data arrays, uses raw format
    template<class T>
    class NakedRawDataArrayWriter : public DataArrayWriter<T>
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
                              int nitems)
        : s(theStream)
      {
        s.write((unsigned int)(ncomps*nitems*sizeof(T)));
      }

      //! write one data element to output stream
      void write (T data)
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
      unsigned bytecount;
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
        : type(type_), stream(stream_), bytecount(0), phase(main)
      { }

      //! signal start of the appeneded section
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
        case binary :         return false;
        case binaryappended : return true;
        }
        DUNE_THROW(IOError, "Dune::VTK::DataArrayWriter: unsupported "
                   "OutputType " << type);
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
       *
       * The should never be more than one DataArrayWriter on the same stream
       * around.  The returned object should be freed with delete.
       */
      template<typename T>
      DataArrayWriter<T>* make(const std::string& name, unsigned ncomps,
                               unsigned nitems, const Indent& indent) {
        switch(phase) {
        case main :
          switch(type) {
          case ascii :
            return new AsciiDataArrayWriter<T>(stream, name, ncomps, indent);
          case binary :
            return new BinaryDataArrayWriter<T>(stream, name, ncomps, nitems,
                                                indent);
          case binaryappended :
            return new BinaryAppendedDataArrayWriter<T>(stream, name, ncomps,
                                                        bytecount, indent);
          }
          break;
        case appended :
          switch(type) {
          case ascii :
          case binary :
            break; // invlid in appended mode
          case binaryappended :
            return new NakedRawDataArrayWriter<T>(stream, ncomps, nitems);
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
