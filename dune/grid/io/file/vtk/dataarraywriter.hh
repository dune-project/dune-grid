// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH

#include <ostream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/vtk/b64enc.hh>
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

  //! base class for data array writers
  /**
   * \tparam T Type of the data elements to write
   *
   * The usage pattern for this class is as follows:
   * \code
     {
     // create an automatically cleaned up writer, giving it the stream to write
     // to, the name the array should have in the VTK file, the number of
     // components (dimension) of the field, and the number of entries in the
     // array overall (number of components times number of cells/number of
     // vertices).  The type of the writer (ascii, binary, binaryappended) is
     // determined from the type argument.
     std::auto_ptr<VTKDataArrayWriter<T> > p
     (makeVTKDataArrayWriter<T>(type, stream, func.name(), func.ncomps(),
                               array_size));
     // Iterate over Elements
     // For point data iterate over the vertices instead (not shown here)
     for (CellIterator i=cellBegin(); i!=cellEnd(); ++i) {
     // iterate over the components
     for (int j=0; j < func.ncomps(); j++)
      p->write(func.evaluate(j, *i, i.position()));
     // Fill 2D data to 3D; makes paraview happier
     if(func.ncomps()==2)
      p->write(0.0);
     }
     }
   * \endcode
   *
   * In the constructor, the actual writer implementation will prepare the
   * stream for writing the data.  This usually means writing a "<DataArray>"
   * header to the stream.  In the write() method it will write a data element
   * to the stream (this is not true for binaryappended mode however: in this
   * mode write just counts the number of bytes that would be written, the
   * actual writing has to happen later independent of the writer).  Finally,
   * in the destructor, the stream is put back in a sane state.  That usually
   * means writing something line "</DataArray>".

   * This is an abstract base class; for an actual implementation look at
   * VTKAsciiDataArrayWriter, VTKBinaryDataArrayWriter, or
   * VTKBinaryAppendedDataArrayWriter.
   */
  template<class T>
  class VTKDataArrayWriter
  {
  public:
    //! write one data element
    virtual void write (T data) = 0;
    //! virtual destructor
    virtual ~VTKDataArrayWriter () {}
  };

  //! a streaming writer for data array tags, uses ASCII inline format
  template<class T>
  class VTKAsciiDataArrayWriter : public VTKDataArrayWriter<T>
  {
  public:
    //! make a new data array writer
    /**
     * \param theStream Stream to write to.
     * \param name      Name of array to write.
     * \param ncomps    Number of components of the array.
     */
    VTKAsciiDataArrayWriter (std::ostream& theStream, std::string name,
                             int ncomps)
      : s(theStream), counter(0), numPerLine(12)
    {
      VTKTypeNameTraits<T> tn;
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
      typedef typename VTKTypeNameTraits<T>::PrintType PT;
      s << (PT) data << " ";
      counter++;
      if (counter%numPerLine==0) s << "\n";
    }

    //! finish output; writes end tag
    ~VTKAsciiDataArrayWriter ()
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
  class VTKBinaryDataArrayWriter : public VTKDataArrayWriter<T>
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
    VTKBinaryDataArrayWriter (std::ostream& theStream, std::string name,
                              int ncomps, int nitems)
      : s(theStream)
    {
      VTKTypeNameTraits<T> tn;
      ncomps = (ncomps>1 ? 3 : 1);
      s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
      // vtk file format: a vector data always should have 3 comps (with 3rd
      // comp = 0 in 2D case)
      if (ncomps>3)
        DUNE_THROW(IOError, "VTKWriter does not support more than 3 "
                   "components");
      s << "NumberOfComponents=\"" << ncomps << "\" ";
      s << "format=\"binary\">\n";
      // reset chunk
      chunk.txt.read(0,0);
      // store size
      unsigned long int size = ncomps*nitems*sizeof(T);
      b64enc(size);
      flush();
    }

    //! write one data element to output stream
    void write (T data)
    {
      b64enc(data);
    }

    //! finish output; writes end tag
    ~VTKBinaryDataArrayWriter ()
    {
      flush();
      s << "\n</DataArray>\n";
      s.flush();
    }

  private:
    template <class X>
    void b64enc(X & data)
    {
      char* p = reinterpret_cast<char*>(&data);
      for (size_t len = sizeof(X); len > 0; len--,p++)
      {
        chunk.txt.put(*p);
        if (chunk.txt.size == 3)
        {
          chunk.data.write(obuf);
          s.write(obuf,4);
        }
      }
    }

    void flush()
    {
      if (chunk.txt.size > 0)
      {
        chunk.data.write(obuf);
        s.write(obuf,4);
      }
    }

    std::ostream& s;
    b64chunk chunk;
    char obuf[4];
  };

  //! a streaming writer for data array tags, uses binary appended format
  template<class T>
  class VTKBinaryAppendedDataArrayWriter : public VTKDataArrayWriter<T>
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
    VTKBinaryAppendedDataArrayWriter (std::ostream& theStream,
                                      std::string name, int ncomps,
                                      unsigned int& bc)
      : s(theStream),bytecount(bc)
    {
      VTKTypeNameTraits<T> tn;
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

  /** @brief Make a VTKDataArrayWriter with new
   *
   * @param outputtype  Which incarnation of VTKDataArrayWriter to return.
   * @param s           The stream to write to
   * @param name        The name of the vtk array
   * @param components  The number of components of the vector
   * @param nitems      Number of vectors in the array (i.e. number of
   *                    cells/vertices)
   * @param bc          Byte count variable: this is incremented by one for
   *                    each byte which has to beritte to the appended data
   *                    section later.  This will actually only happen for
   *                    outputtype=binaryappended.
   */
  template<class T> VTKDataArrayWriter<T> *
  makeVTKDataArrayWriter(VTKOptions::OutputType outputtype, std::ostream &s,
                         const char *name, unsigned int components,
                         unsigned int nitems, unsigned int& bc)
  {
    switch(outputtype) {
    case VTKOptions::ascii :
      return new VTKAsciiDataArrayWriter<T>(s, name, components);
    case VTKOptions::binary :
      return new VTKBinaryDataArrayWriter<T>(s, name, components, nitems);
    case VTKOptions::binaryappended :
      return new VTKBinaryAppendedDataArrayWriter<T>(s, name, components, bc);
    }
    DUNE_THROW(IOError, "VTKDataArrayWriter: unsupported OutputType " <<
               outputtype);
  }

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_DATAARRAYWRITER_HH
