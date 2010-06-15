// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH

#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>

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
       */
      AsciiDataArrayWriter(std::ostream& theStream, std::string name,
                           int ncomps)
        : s(theStream), counter(0), numPerLine(12)
      {
        TypeName<T> tn;
        s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
        // vtk file format: a vector data always should have 3 comps (with 3rd
        // comp = 0 in 2D case)
        if (ncomps>3)
          DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                     "components");
        s << "NumberOfComponents=\"" << (ncomps>1 ? 3 : 1) << "\" ";
        s << "format=\"ascii\">\n";
      }

      //! write one data element to output stream
      void write (T data)
      {
        typedef typename PrintType<T>::Type PT;
        s << (PT) data << " ";
        counter++;
        if (counter%numPerLine==0) s << "\n";
      }

      //! finish output; writes end tag
      ~AsciiDataArrayWriter ()
      {
        if (counter%numPerLine!=0) s << std::endl;
        s << "</DataArray>\n";
      }

    private:
      std::ostream& s;
      int counter;
      int numPerLine;
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
       */
      BinaryDataArrayWriter(std::ostream& theStream, std::string name,
                            int ncomps, int nitems)
        : s(theStream), b64(theStream)
      {
        TypeName<T> tn;
        ncomps = (ncomps>1 ? 3 : 1);
        s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
        // vtk file format: a vector data always should have 3 comps (with 3rd
        // comp = 0 in 2D case)
        if (ncomps>3)
          DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                     "components");
        s << "NumberOfComponents=\"" << ncomps << "\" ";
        s << "format=\"binary\">\n";

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
        s << "\n</DataArray>\n";
        s.flush();
      }

    private:
      std::ostream& s;
      Base64Stream b64;
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
       */
      BinaryAppendedDataArrayWriter(std::ostream& theStream, std::string name,
                                    int ncomps, unsigned int& bc)
        : s(theStream),bytecount(bc)
      {
        TypeName<T> tn;
        s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
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

    //! a factory for DataArrayWriters
    /**
     * Some types of DataArrayWriters need to communicate data sauch as byte
     * counts between different instances which write to the same stream.
     * This factory will manage that data.
     */
    class DataArrayWriterFactory {
      OutputType type;
      std::ostream& stream;
      unsigned bytecount;

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
        : type(type_), stream(stream_), bytecount(0)
      { }

      //! create a DataArrayWriter
      /**
       * \tparam T Type of the data to write.
       *
       * \param name   Name of the array to write.
       * \param ncomps Number of components of the vectors in to array.
       * \param nitems Number of vectors in the array.
       *
       * The should never be more than one DataArrayWriter on the same stream
       * around.  The returned object should be freed with delete.
       */
      template<typename T>
      DataArrayWriter<T>* make(const std::string& name, unsigned ncomps,
                               unsigned nitems) {
        switch(type) {
        case ascii :
          return new AsciiDataArrayWriter<T>(stream, name, ncomps);
        case binary :
          return new BinaryDataArrayWriter<T>(stream, name, ncomps, nitems);
        case binaryappended :
          return new BinaryAppendedDataArrayWriter<T>(stream, name, ncomps,
                                                      bytecount);
        }
        DUNE_THROW(IOError, "Dune::VTK::DataArrayWriter: unsupported "
                   "OutputType " << type);
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
