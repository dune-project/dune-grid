// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_STREAMS_HH
#define DUNE_GRID_IO_FILE_VTK_STREAMS_HH

#include <ostream>

#include <dune/grid/io/file/vtk/b64enc.hh>

namespace Dune {

  //! class to base64 encode a stream of data
  class Base64Stream {
    std::ostream& s;
    b64chunk chunk;
    char obuf[4];

  public:
    //! Construct a Base64Stream
    /**
     * \param s_ The stream the resulting base64-encoded text will be written
     *           to.
     */
    Base64Stream(std::ostream& s_)
      : s(s_)
    {
      // reset chunk
      chunk.reset();
    }

    //! encode a data item
    /**
     * The result will be written to the stream, eventually.  This method may
     * be called multiple times in a row.  After this method has been called,
     * no one else may write to the underlying stream until flush() has been
     * called or this writer object has been destroyed.
     */
    template <class X>
    void write(X & data)
    {
      char* p = reinterpret_cast<char*>(&data);
      for (size_t len = sizeof(X); len > 0; len--,p++)
      {
        chunk.put(*p);
        if (chunk.size == 3)
        {
          chunk.write(obuf);
          s.write(obuf,4);
        }
      }
    }

    //! flush the current unwritten data to the stream.
    /**
     * If the size of the received input is not a multiple of three bytes, an
     * end-marker will be written.
     *
     * Calling this function a second time without calling b64enc() or calling
     * it right after construction has no effect.
     */
    void flush()
    {
      if (chunk.size > 0)
      {
        chunk.write(obuf);
        s.write(obuf,4);
      }
    }

    //! destroy the object
    /**
     * Calls flush()
     */
    ~Base64Stream() {
      flush();
    }
  };

  //! write out data in binary
  class RawStream
  {
  public:
    //! make a new stream
    inline RawStream (std::ostream& theStream)
      : s(theStream)
    {}

    //! write data to stream
    template<class T>
    void write (T data)
    {
      char* p = reinterpret_cast<char*>(&data);
      s.write(p,sizeof(T));
    }
  private:
    std::ostream& s;
  };

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_STREAMS_HH
