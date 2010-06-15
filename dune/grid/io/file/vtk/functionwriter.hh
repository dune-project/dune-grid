// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH

#include <cstddef>
#include <string>
#include <typeinfo>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
#include <dune/grid/io/file/vtk/pvtuwriter.hh>
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

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) = 0;
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

    //////////////////////////////////////////////////////////////////////
    //
    //  A Generic Function writer for VTKFunctions
    //

    //! Base class for function writers
    template<typename Func>
    class VTKFunctionWriter
      : public FunctionWriterBase<typename Func::Entity>
    {
      typedef FunctionWriterBase<typename Func::Entity> Base;
      shared_ptr<const Func> func;
      shared_ptr<DataArrayWriter<float> > arraywriter;

    public:
      VTKFunctionWriter(const shared_ptr<const Func>& func_)
        : func(func_)
      { }

      //! return name
      virtual std::string name() const { return func->name(); }

      //! return number of components of the vector
      virtual unsigned ncomps() const {
        if(func->ncomps() == 2) return 3;
        else return func->ncomps();
      }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<float>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<float>(name(), ncomps(),
                                                        nitems));
        return !arraywriter->writeIsNoop();
      }

      //! write at the given position
      virtual void write(const typename Base::Cell& cell,
                         const typename Base::Domain& xl) {
        for(int d = 0; d < func->ncomps(); ++d)
          arraywriter->write(func->evaluate(d, cell, xl));
        for(unsigned d = func->ncomps(); d < ncomps(); ++d)
          arraywriter->write(0);
      }

      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  Writers for the grid information
    //

    //! writer for the Coordinates array
    template<typename Cell>
    class CoordinatesWriter
      : public FunctionWriterBase<Cell>
    {
      typedef FunctionWriterBase<Cell> Base;

      shared_ptr<DataArrayWriter<float> > arraywriter;

    public:

      //! return name
      virtual std::string name() const { return "Coordinates"; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return 3; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<float>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<float>(name(), ncomps(),
                                                        nitems));
        return !arraywriter->writeIsNoop();
      }
      //! write at the given position
      virtual void write(const typename Base::Cell& cell,
                         const typename Base::Domain& xl) {
        FieldVector<typename Base::Cell::ctype, Base::Cell::dimensionworld> xg
          = cell.geometry().global(xl);
        for(unsigned d = 0; d < 3 && d < Base::Cell::dimensionworld; ++d)
          arraywriter->write(xg[d]);
        for(unsigned d = Base::Cell::dimensionworld; d < 3; ++d)
          arraywriter->write(0);
      };
      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

    //! writer for the connectivity array in conforming mode
    template<typename IteratorFactory>
    class ConformingConnectivityWriter
      : public FunctionWriterBase<typename IteratorFactory::Cell>
    {
      typedef FunctionWriterBase<typename IteratorFactory::Cell> Base;
      static const unsigned mydim = Base::Cell::mydimension;

      const IteratorFactory& factory;
      shared_ptr<DataArrayWriter<unsigned> > arraywriter;
      std::vector<unsigned> pointIndices;

    public:
      //! create a writer with the given iteratorfactory
      ConformingConnectivityWriter(const IteratorFactory& factory_)
        : factory(factory_)
      { }

      //! return name
      virtual std::string name() const { return "connectivity"; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return 1; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<unsigned>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<unsigned>(name(), ncomps(),
                                                           nitems));
        if(arraywriter->writeIsNoop())
          return false;

        //! write is meaningful, we need to build the data
        pointIndices.resize(factory.indexSet().size(mydim));
        const typename IteratorFactory::PointIterator& pend =
          factory.endPoints();
        typename IteratorFactory::PointIterator pit = factory.beginPoints();
        unsigned counter = 0;
        while(pit != pend) {
          pointIndices[factory.indexSet().subIndex
                         (pit->cell(), pit->duneIndex(), mydim)] = counter;
          ++counter;
          ++pit;
        }
        return true;
      }
      //! write at the given corner
      virtual void write(const typename Base::Cell& cell, unsigned cornerIndex)
      {
        // if pointIndices is empty, we're in writeIsNoop mode
        if(pointIndices.size() == 0)
          return;
        arraywriter->write(pointIndices[factory.indexSet().subIndex
                                          (cell, cornerIndex, mydim)]);
      }
      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
        pointIndices.clear();
      }
    };

    //! writer for the connectivity array in nonconforming mode
    template<typename Cell>
    class NonConformingConnectivityWriter
      : public FunctionWriterBase<Cell>
    {
      shared_ptr<DataArrayWriter<unsigned> > arraywriter;
      unsigned counter;

    public:
      //! return name
      virtual std::string name() const { return "connectivity"; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return 1; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<unsigned>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<unsigned>(name(), ncomps(),
                                                           nitems));
        counter = 0;
        return !arraywriter->writeIsNoop();
      }
      //! write at the given corner
      virtual void write(const Cell& cell, unsigned cornerIndex)
      {
        arraywriter->write(counter);
        ++counter;
      }
      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

    //! writer for the offsets array
    template<typename Cell>
    class OffsetsWriter
      : public FunctionWriterBase<Cell>
    {
      typedef FunctionWriterBase<Cell> Base;

      shared_ptr<DataArrayWriter<unsigned> > arraywriter;
      unsigned offset;

    public:
      //! return name
      virtual std::string name() const { return "offsets"; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return 1; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<unsigned>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<unsigned>(name(), ncomps(),
                                                           nitems));
        offset = 0;
        return !arraywriter->writeIsNoop();
      }
      //! write at the given position
      virtual void write(const Cell& cell, const typename Base::Domain&) {
        offset += cell.geometry().corners();
        arraywriter->write(offset);
      };
      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

    //! writer for the types array
    template<typename Cell>
    class TypesWriter
      : public FunctionWriterBase<Cell>
    {
      typedef FunctionWriterBase<Cell> Base;

      shared_ptr<DataArrayWriter<unsigned char> > arraywriter;

    public:
      //! return name
      virtual std::string name() const { return "types"; }

      //! return number of components of the vector
      virtual unsigned ncomps() const { return 1; }

      //! add this field to the given parallel writer
      virtual void addArray(PVTUWriter& writer) {
        writer.addArray<unsigned char>(name(), ncomps());
      }

      //! start writing with the given writer
      virtual bool beginWrite(VTUWriter& writer, std::size_t nitems) {
        arraywriter.reset(writer.makeArrayWriter<unsigned char>
                            ( name(), ncomps(), nitems));
        return !arraywriter->writeIsNoop();
      }
      //! write at the given position
      virtual void write(const Cell& cell, const typename Base::Domain&) {
        arraywriter->write(geometryType(cell.type()));
      };
      //! signal end of writing
      virtual void endWrite() {
        arraywriter.reset();
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_FUNCTIONWRITER_HH
