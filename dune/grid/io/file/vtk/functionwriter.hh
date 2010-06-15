// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH

#include <cstddef>
#include <string>
#include <typeinfo>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! Base class for function writers
    template<typename Cell_>
    class FunctionWriterBase {
      typedef typename Cell_::ctype DF;
      static const unsigned mydim = Cell_::mydimension;
      typedef GenericReferenceElements<DF, mydim> Refelems;

    public:
      typedef FieldVector<DF, mydim> Domain;
      typedef Cell_ Cell;

      //! return name
      virtual std::string name() const = 0;

      //! return number of components of the vector
      virtual unsigned ncomps() const = 0;

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) = 0;
      //! write at the given position
      /**
       * This is the default dummy implementation.  This method is not
       * abstract so derived classes don't have to override it if they don't
       * need it.
       */
      virtual void write(const Cell& cell, const Domain& xl) {
        DUNE_THROW(NotImplemented, "FunctionWriterBase::write(const Cell&, "
                   "const Domain&): Either the derived class " <<
                   typeid(*this).name() << " failed to implement this method "
                   "or this method is not meant to be called on the derived "
                   "class and was called in error.");
      };
      //! write at the given corner
      /**
       * This default method forwards the writing to write(const Cell&, const
       * Domain&).
       */
      virtual void write(const Cell& cell, unsigned cornerIndex) {
        write(cell,
              Refelems::general(cell.type()).position(cornerIndex, mydim));
      }
      //! signal end of writing
      virtual void endWrite() = 0;
      //! destructor
      virtual ~FunctionWriterBase() {};
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH
