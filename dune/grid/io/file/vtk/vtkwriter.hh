// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_VTKWRITER_HH
#define DUNE_VTKWRITER_HH

#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <type_traits>
#include <vector>
#include <list>
#include <map>

#include <dune/common/visibility.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/path.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/pvtuwriter.hh>
#include <dune/grid/io/file/vtk/streams.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Provides file i/o for the visualization toolkit
 */

/**
   @todo put vtk io intro here ...

   details and examples regarding the VTK file format can be found here:

   http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats
*/

namespace Dune
{

  namespace Impl
  {
    // Check whether type F has a method 'bind'  (see the dune-functions interface)
    template< class F, class E, class = void >
    struct IsBindable
      : std::false_type
    {};

    template< class F, class E >
    struct IsBindable< F, E, std::void_t< decltype( std::declval< F & >().bind( std::declval< const E & >() ) ),
                                          decltype( std::declval< F & >().unbind() ) > >
      : std::true_type
    {};

    // Check whether localFunction(F) can be called  (see the dune-functions interface)
    template< class F, class = void >
    struct HasLocalFunction
      : std::false_type
    {};

    template< class F >
    struct HasLocalFunction< F, std::void_t< decltype( localFunction( std::declval< F& >() ) ) > >
      : std::true_type
    {};

  } // namespace Impl

  // Forward-declaration here, so the class can be friend of VTKWriter
  template <class GridView>
  class VTKSequenceWriterBase;
  template <class GridView>
  class VTKSequenceWriter;

  /**
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class VTKWriter {

    // VTKSequenceWriterBase needs getSerialPieceName
    // and getParallelHeaderName
    friend class VTKSequenceWriterBase<GridView>;
    // VTKSequenceWriter needs the grid view, to get the MPI size and rank
    friend class VTKSequenceWriter<GridView>;

    // extract types
    typedef typename GridView::Grid Grid;
    typedef typename GridView::ctype DT;
    constexpr static int n = GridView::dimension;
    constexpr static int w = GridView::dimensionworld;

    typedef typename GridView::template Codim< 0 >::Entity Cell;
    typedef typename GridView::template Codim< n >::Entity Vertex;
    typedef Cell Entity;

    typedef typename GridView::IndexSet IndexSet;

    static const PartitionIteratorType VTK_Partition = InteriorBorder_Partition;
    //static const PartitionIteratorType VTK_Partition = All_Partition;

    typedef typename GridView::template Codim< 0 >
    ::template Partition< VTK_Partition >::Iterator
    GridCellIterator;
    typedef typename GridView::template Codim< n >
    ::template Partition< VTK_Partition >::Iterator
    GridVertexIterator;

    typedef typename GridCellIterator::Reference EntityReference;

    typedef typename GridView::template Codim< 0 >
    ::Entity::Geometry::LocalCoordinate Coordinate;

    typedef MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    // return true if entity should be skipped in Vertex and Corner iterator
    static bool skipEntity( const PartitionType entityType )
    {
      switch( VTK_Partition )
      {
        // for All_Partition no entity has to be skipped
        case All_Partition:             return false;
        case InteriorBorder_Partition:  return ( entityType != InteriorEntity );
        default: DUNE_THROW(NotImplemented,"Add check for this partition type");
      }
      return false ;
    }

  public:

    typedef Dune::VTKFunction< GridView > VTKFunction;

  protected:

    //! Type erasure wrapper for VTK data sets
    /**
     * This wrapper has value semantics
     */
    class VTKLocalFunction
    {

    public:

      typedef VTK::DataArrayWriter Writer;

      //! Base class for polymorphic container of underlying data set
      struct FunctionWrapperBase
      {

        //! Bind data set to grid entity - must be called before evaluating (i.e. calling write())
        virtual void bind(const Entity& e) = 0;

        //! Unbind data set from current grid entity - mostly here for performance and symmetry reasons
        virtual void unbind() = 0;

        //! Evaluate data set at local position pos inside the current entity and write result to w.
        /**
         * The function must write count scalar values as determined by the VTK::FieldInfo.
         */
        virtual void write(const Coordinate& pos, Writer& w, std::size_t count) const = 0;

        virtual ~FunctionWrapperBase()
        {}

      };

      //! Type erasure implementation for functions conforming to the dune-functions LocalFunction interface
      // DUNE_PRIVATE since _f has less visibility
      template<typename F>
      struct DUNE_PRIVATE FunctionWrapper
        : public FunctionWrapperBase
      {
        using Function = typename std::decay<F>::type;

        template<typename F_>
        FunctionWrapper(F_&& f)
          : _f(std::forward<F_>(f))
        {}

        virtual void bind(const Entity& e)
        {
          _f.bind(e);
        }

        virtual void unbind()
        {
          _f.unbind();
        }

        virtual void write(const Coordinate& pos, Writer& w, std::size_t count) const
        {
          auto r = _f(pos);
          // we need to do different things here depending on whether r supports indexing into it or not.
          do_write(w,r,count,IsIndexable<decltype(r)>());
        }

      private:

        template<typename R>
        void do_write(Writer& w, const R& r, std::size_t count, std::true_type) const
        {
          for (std::size_t i = 0; i < count; ++i)
            w.write(r[i]);
        }

        template<typename R>
        void do_write(Writer& w, const R& r, std::size_t count, std::false_type) const
        {
          assert(count == 1);
          w.write(r);
        }

        Function _f;
      };

      //! Type erasure implementation for C++ functions, i.e., functions that can be evaluated in global coordinates
      template<typename F>
      struct GlobalFunctionWrapper
        : public FunctionWrapperBase
      {
        using Function = typename std::decay<F>::type;

        template<typename F_>
        GlobalFunctionWrapper(F_&& f)
          : _f(std::forward<F_>(f))
          , element_(nullptr)
        {}

        virtual void bind(const Entity& e)
        {
          element_ = &e;
        }

        virtual void unbind()
        {
          element_ = nullptr;
        }

        virtual void write(const Coordinate& pos, Writer& w, std::size_t count) const
        {
          auto globalPos = element_->geometry().global(pos);
          auto r = _f(globalPos);
          if constexpr (IsIndexable<decltype(r)>()) {
            for (std::size_t i = 0; i < count; ++i)
                w.write(r[i]);
          }
          else {
            assert(count == 1);
            w.write(r);
          }
        }
      private:
        Function _f;
        const Entity* element_;
      };

      //! Type erasure implementation for legacy VTKFunctions.
      struct VTKFunctionWrapper
        : public FunctionWrapperBase
      {
        VTKFunctionWrapper(const std::shared_ptr< const VTKFunction >& f)
          : _f(f)
          , _entity(nullptr)
        {}

        virtual void bind(const Entity& e)
        {
          _entity = &e;
        }

        virtual void unbind()
        {
          _entity = nullptr;
        }

        virtual void write(const Coordinate& pos, Writer& w, std::size_t count) const
        {
          for (std::size_t i = 0; i < count; ++i)
            w.write(_f->evaluate(i,*_entity,pos));
        }

      private:

        std::shared_ptr< const VTKFunction > _f;
        const Entity* _entity;

      };

      //! Construct a VTKLocalFunction for a dune-functions style LocalFunction
      template<typename F, std::enable_if_t<Impl::IsBindable<F, Entity>::value, int> = 0>
      VTKLocalFunction(F&& f, VTK::FieldInfo fieldInfo)
        : _f(std::make_unique<FunctionWrapper<F> >(std::forward<F>(f)))
        , _fieldInfo(fieldInfo)
      {}

      //! Construct a VTKLocalFunction for a dune-functions GridViewFunction
      // That is, a function that you can create a LocalFunction for, and evaluate that in element coordinates
      template<typename F, std::enable_if_t<not Impl::IsBindable<F, Entity>::value && Impl::HasLocalFunction<F>::value, int> = 0>
      VTKLocalFunction(F&& f, VTK::FieldInfo fieldInfo)
        : _f(std::make_unique< FunctionWrapper<
          typename std::decay<decltype(localFunction(std::forward<F>(f)))>::type
          > >(localFunction(std::forward<F>(f))))
        , _fieldInfo(fieldInfo)
      {}

      //! Construct a VTKLocalFunction for a C++ (global) function
      // That is, a function that can be evaluated in global coordinates of the domain
      template<typename F, std::enable_if_t<not Impl::IsBindable<F, Entity>::value && not Impl::HasLocalFunction<F>::value, int> = 0>
      VTKLocalFunction(F&& f, VTK::FieldInfo fieldInfo)
        : _f(std::make_unique< GlobalFunctionWrapper<F> >(std::forward<F>(f)))
        , _fieldInfo(fieldInfo)
      {}

      //! Construct a VTKLocalFunction for a legacy VTKFunction
      explicit VTKLocalFunction (const std::shared_ptr< const VTKFunction >& vtkFunctionPtr)
        : _f(std::make_unique<VTKFunctionWrapper>(vtkFunctionPtr))
        , _fieldInfo(
          vtkFunctionPtr->name(),
          (vtkFunctionPtr->ncomps() == 2 || vtkFunctionPtr->ncomps() == 3)  ? VTK::FieldInfo::Type::vector : VTK::FieldInfo::Type::scalar,
          vtkFunctionPtr->ncomps(),
          vtkFunctionPtr->precision()
          )
      {}

      //! Returns the name of the data set
      std::string name() const
      {
        return fieldInfo().name();
      }

      //! Returns the VTK::FieldInfo for the data set
      const VTK::FieldInfo& fieldInfo() const
      {
        return _fieldInfo;
      }

      //! Bind the data set to grid entity e.
      void bind(const Entity& e) const
      {
        _f->bind(e);
      }

      //! Unbind the data set from the currently bound entity.
      void unbind() const
      {
        _f->unbind();
      }

      //! Write the value of the data set at local coordinate pos to the writer w.
      void write(const Coordinate& pos, Writer& w) const
      {
        _f->write(pos,w,fieldInfo().size());
      }

      std::shared_ptr<FunctionWrapperBase> _f;
      VTK::FieldInfo _fieldInfo;

    };

    typedef typename std::list<VTKLocalFunction>::const_iterator FunctionIterator;

    //! Iterator over the grids elements
    /**
     * This class iterates over the gridview's elements.  It is the same as
     * the gridview's Codim<0>::Iterator for the InteriorBorder_Partition,
     * except that it add a position() method.
     */
    class CellIterator : public GridCellIterator
    {
    public:
      //! construct a CellIterator from the gridview's Iterator.
      CellIterator(const GridCellIterator & x) : GridCellIterator(x) {}
      //! get the position of the center of the element, in element-local
      //! coordinates
      const FieldVector<DT,n> position() const
      {
        return ReferenceElements<DT,n>::general((*this)->type()).position(0,0);
      }
    };

    CellIterator cellBegin() const
    {
      return gridView_.template begin< 0, VTK_Partition >();
    }

    CellIterator cellEnd() const
    {
      return gridView_.template end< 0, VTK_Partition >();
    }

    //! Iterate over the grid's vertices
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  If the data mode dm is nonconforming, each vertex is visited
     * once for each element where it is a corner (similar to CornerIterator).
     * If dm is conforming each vertex is visited only once globally, for the
     * first element where it is a corner.  Contrary to CornerIterator, visit
     * the corners of a given element in Dune-ordering.
     *
     * Dereferencing the iterator yields the current entity, and the index of
     * the current corner within that entity is returned by the iterators
     * localindex() method.  Another useful method on the iterator itself is
     * position() which returns the element-local position of the current
     * corner.
     */
    class VertexIterator :
      public ForwardIteratorFacade<VertexIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in Dune-numbering, in contrast to CornerIterator.
      int cornerIndexDune;
      const VertexMapper & vertexmapper;
      std::vector<bool> visited;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order (VertexIterator)
      int offset;

      // hide operator ->
      void operator->();
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++cornerIndexDune;
        const int numCorners = git->subEntities(n);
        if( cornerIndexDune == numCorners )
        {
          offset += numCorners;
          cornerIndexDune = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
    public:
      VertexIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const VTK::DataMode & dm,
                     const VertexMapper & vm) :
        git(x), gend(end), datamode(dm), cornerIndexDune(0),
        vertexmapper(vm), visited(vm.size(), false),
        offset(0)
      {
        if (datamode == VTK::conforming && git != gend)
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
      }
      void increment ()
      {
        switch (datamode)
        {
        case VTK::conforming :
          while(visited[vertexmapper.subIndex(*git,cornerIndexDune,n)])
          {
            basicIncrement();
            if (git == gend) return;
          }
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
          break;
        case VTK::nonconforming :
          basicIncrement();
          break;
        }
      }
      bool equals (const VertexIterator & cit) const
      {
        return git == cit.git
               && cornerIndexDune == cit.cornerIndexDune
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! index of vertex within the entity, in Dune-numbering
      int localindex () const
      {
        return cornerIndexDune;
      }
      //! position of vertex inside the entity
      FieldVector<DT,n> position () const
      {
        return referenceElement<DT,n>(git->type())
          .position(cornerIndexDune,n);
      }
    };

    VertexIterator vertexBegin () const
    {
      return VertexIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    VertexIterator vertexEnd () const
    {
      return VertexIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    //! Iterate over the elements' corners
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  Each vertex in the grid can be a corner in multiple elements,
     * and is visited once for each element it is associated with.  This class
     * differs from VertexIterator in that it visits the corners of a given
     * element in VTK-ordering, and that it always visits a given vertex once
     * for each element where that vertex is a corner in, independent of the
     * data mode dm.
     *
     * Dereferencing the iterator yields the current entity.  Another useful
     * method on the iterator itself is id(), which returns the number of the
     * current corners associated vertex, in the numbering given by the
     * iteration order of VertexIterator.
     */
    class CornerIterator :
      public ForwardIteratorFacade<CornerIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in VTK-numbering, in contrast to VertexIterator.
      int cornerIndexVTK;
      const VertexMapper & vertexmapper;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order of VertexIterator (*not*
      // CornerIterator)
      const std::vector<int> & number;
      // holds the number of corners of all the elements we have seen so far,
      // excluding the current element
      int offset;

      // hide operator ->
      void operator->();
    public:
      CornerIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const VTK::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), cornerIndexVTK(0),
        vertexmapper(vm),
        number(num), offset(0) {}
      void increment ()
      {
        if( git == gend )
          return;
        ++cornerIndexVTK;
        const int numCorners = git->subEntities(n);
        if( cornerIndexVTK == numCorners )
        {
          offset += numCorners;
          cornerIndexVTK = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
      bool equals (const CornerIterator & cit) const
      {
        return git == cit.git
               && cornerIndexVTK == cit.cornerIndexVTK
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! Process-local consecutive zero-starting vertex id
      /**
       * This method returns the number of this corners associated vertex, in
       * the numbering given by the iteration order of VertexIterator.
       */
      int id () const
      {
        switch (datamode)
        {
        case VTK::conforming :
          return
            number[vertexmapper.subIndex(*git,VTK::renumber(*git,cornerIndexVTK),
                                    n)];
        case VTK::nonconforming :
          return offset + VTK::renumber(*git,cornerIndexVTK);
        default :
          DUNE_THROW(IOError,"VTKWriter: unsupported DataMode" << datamode);
        }
      }
    };

    CornerIterator cornerBegin () const
    {
      return CornerIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

    CornerIterator cornerEnd () const
    {
      return CornerIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

  public:
    /**
     * @brief Construct a VTKWriter working on a specific GridView.
     *
     *
     * @param gridView The gridView the grid functions live on. (E. g. a LevelGridView.)
     * @param dm The data mode.
     * @param coordPrecision the precision with which to write out the coordinates
     */
    explicit VTKWriter ( const GridView &gridView,
                         VTK::DataMode dm = VTK::conforming,
                         VTK::Precision coordPrecision = VTK::Precision::float32)
      : gridView_( gridView ),
        datamode( dm ),
        coordPrec (coordPrecision),
        polyhedralCellsPresent_( checkForPolyhedralCells() )
    { }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p std::shared_ptr to the function to visualize
     */
    void addCellData (const std::shared_ptr< const VTKFunction > & p)
    {
      celldata.push_back(VTKLocalFunction(p));
    }

    /**
     * @brief Add a function by sampling it on the element centers
     *
     * \param f The function to be written to the file
     *
     * The object f can be one of several things.  Depending on what it is exactly,
     * its object life-time is influenced in different ways:
     * - If f has the method bind(), then a copy of f is stored, regardless of whether f is an l- or r-value.
     * - If f can be localized by calling localFunction(f), then a copy of localFunction(f) is stored,
     *   but f is never stored, regardless of whether f is an l- or r-value.
     * - If f supports neither bind() or localFunction(), then a copy of f is stored,
     *   regardless if f is an l- or r-value.
     *
     * The previous paragraph actually refers to parts of the dune-functions interface,
     * and you may want to read up on that if you want to write functions to VTK.
     *
     * \deprecated f may also be a VTKFunction object, but you are strongly discouraged
     *   from using VTKFunctions.
     */
    template<typename F>
    void addCellData(F&& f, VTK::FieldInfo vtkFieldInfo)
    {
      celldata.push_back(VTKLocalFunction(std::forward<F>(f),vtkFieldInfo));
    }

    /**
     * @brief Add a grid function (represented by container) that lives on the cells of
     * the grid to the visualization.
     *
     * The container has to have random access via operator[] (e.g. std::vector). The
     * value of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the index (corresponding
     * to the index from the MGMC mapper on the grid view) of the element.
     * For vector valued data all components for an element are assumed to
     * be consecutive.
     *
     * @param v The container with the values of the grid function for each cell.
     * @param name A name to identify the grid function.
     * @param ncomps Number of components (default is 1).
     */
    template<class Container>
    void addCellData (const Container& v, const std::string &name, int ncomps = 1,
                      VTK::Precision prec = VTK::Precision::float32)
    {
      typedef P0VTKFunction<GridView, Container> Function;
      for (int c=0; c<ncomps; ++c) {
        std::stringstream compName;
        compName << name;
        if (ncomps>1)
          compName << "[" << c << "]";
        VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c, prec);
        addCellData(std::shared_ptr< const VTKFunction >(p));
      }
    }

    /**
     * @brief Add a grid function that lives on the vertices of the grid to the visualization.
     * @param p std::shared_ptr to the function to visualize
     */
    void addVertexData (const std::shared_ptr< const VTKFunction > & p)
    {
      vertexdata.push_back(VTKLocalFunction(p));
    }

    /**
     * @brief Add a function by sampling it on the grid vertices
     *
     * \param f The function to be written to the file
     *
     * The object f can be one of several things.  Depending on what it is exactly,
     * its object life-time is influenced in different ways:
     * - If f has the method bind(), then a copy of f is stored, regardless of whether f is an l- or r-value.
     * - If f can be localized by calling localFunction(f), then a copy of localFunction(f) is stored,
     *   but f is never stored, regardless of whether f is an l- or r-value.
     * - If f supports neither bind() or localFunction(), then a copy of f is stored,
     *   regardless if f is an l- or r-value.
     *
     * The previous paragraph actually refers to parts of the dune-functions interface,
     * and you may want to read up on that if you want to write functions to VTK.
     *
     * \deprecated f may also be a VTKFunction object, but you are strongly discouraged
     *   from using VTKFunctions.
     */
    template<typename F>
    void addVertexData(F&& f, VTK::FieldInfo vtkFieldInfo)
    {
      vertexdata.push_back(VTKLocalFunction(std::forward<F>(f),vtkFieldInfo));
    }


    /**
     * @brief Add a grid function (represented by container) that lives on the vertices of the
     * grid to the visualization output.
     *
     * The container has to have random access via operator[] (e.g. std::vector). The value
     * of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the index (corresponding
     * to the index from the MGMC mapper on the grid view) of the vertex.
     * For vector valued data all components for a vertex are assumed to
     * be consecutive.
     *
     * @param v The container with the values of the grid function for each vertex.
     * @param name A name to identify the grid function.
     * @param ncomps Number of components (default is 1).
     */
    template<class Container>
    void addVertexData (const Container& v, const std::string &name, int ncomps=1,
                        VTK::Precision prec = VTK::Precision::float32)
    {
      typedef P1VTKFunction<GridView, Container> Function;
      for (int c=0; c<ncomps; ++c) {
        std::stringstream compName;
        compName << name;
        if (ncomps>1)
          compName << "[" << c << "]";
        VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c, prec);
        addVertexData(std::shared_ptr< const VTKFunction >(p));
      }
    }

    //! clear list of registered functions
    void clear ()
    {
      celldata.clear();
      vertexdata.clear();
    }

    //! get the precision with which coordinates are written out
    VTK::Precision coordPrecision() const
    { return coordPrec; }

    //! destructor
    virtual ~VTKWriter ()
    {
      this->clear();
    }

    /** \brief write output (interface might change later)
     *
     *  This method can be used in parallel as well as in serial programs.
     *  For serial runs (commSize=1) it chooses other names without the
     *  "s####-p####-" prefix for the .vtu/.vtp files and omits writing of the
     *  .pvtu/pvtp file however.  For parallel runs (commSize > 1) it is the
     *  same as a call to pwrite() with name and path constructed from
     *  a given filename possibly containing a path, and extendpath="".
     *
     *  \param[in]  name  basic name to write (may not contain a path)
     *  \param[in]  type  type of output (e.g,, ASCII) (optional)
     */
    std::string write ( const std::string &name,
                        VTK::OutputType type = VTK::ascii )
    {
      return write( name, type, gridView_.comm().rank(), gridView_.comm().size() );
    }

    /** \brief write output (interface might change later)
     *
     * "pwrite" means "path write" (i.e. write somewhere else than the current
     * directory).  The "p" does not mean this method has a monopoly on
     * parallel writing, the regular write(const std::string &,
     * VTK::OutputType) method can do that just fine.
     *
     * \param name       Base name of the output files.  This should not
     *                   contain any directory part and not filename
     *                   extensions.  It will be used both for each processes
     *                   piece as well as the parallel collection file.
     * \param path       Directory where to put the parallel collection
     *                   (.pvtu/.pvtp) file.  If it is relative, it is taken
     *                   relative to the current directory.
     * \param extendpath Directory where to put the piece file (.vtu/.vtp) of
     *                   this process.  If it is relative, it is taken
     *                   relative to the directory denoted by path.
     * \param type       How to encode the data in the file.
     *
     * \note Currently, extendpath may not be absolute unless path is
     *       absolute, because that would require the value of the current
     *       directory.
     *
     * \throw NotImplemented Extendpath is absolute but path is relative.
     * \throw IOError        Failed to open a file.
     */
    std::string pwrite ( const std::string & name,  const std::string & path, const std::string & extendpath,
                         VTK::OutputType type = VTK::ascii )
    {
      return pwrite( name, path, extendpath, type, gridView_.comm().rank(), gridView_.comm().size() );
    }

  protected:
    //! return name of a parallel piece file (or header name)
    /**
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting piece name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     * \param commRank Rank of the process to generate a piece name for. if (-1)
     * then the header is created.
     * \param commSize Number of processes writing a parallel vtk output.
     */
    std::string getParallelPieceName(const std::string& name,
                                     const std::string& path,
                                     int commRank, int commSize) const
    {
      std::ostringstream s;
      // write path first
      if(path.size() > 0)
      {
        s << path;
        if(path[path.size()-1] != '/')
          s << '/';
      }

      std::string fileprefix;
      // check if a path was already added to name
      // and if yes find filename without path
      auto pos = name.rfind('/');
      if( pos != std::string::npos )
      {
        // extract filename without path
        fileprefix = name.substr( pos+1 );
        // extract the path and added it before
        // the magic below is added
        std::string newpath = name.substr(0, pos);
        s << newpath;
        if(newpath[name.size()-1] != '/')
          s << '/';
      }
      else
      {
        // if no path was found just copy the name
        fileprefix = name;
      }

      s << 's' << std::setw(4) << std::setfill('0') << commSize << '-';
      const bool writeHeader = commRank < 0;
      if( ! writeHeader )
      {
        s << 'p' << std::setw(4) << std::setfill('0') << commRank << '-';
      }

      s << fileprefix << ".";
      // write p for header files
      if( writeHeader )
        s << "p";
      s << "vt";

      if(GridView::dimension > 1)
        s << "u";
      else
        s << "p";
      return s.str();
    }

    //! return name of a parallel header file
    /**
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting header name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     * \param commSize Number of processes writing a parallel vtk output.
     */
    std::string getParallelHeaderName(const std::string& name,
                                      const std::string& path,
                                      int commSize) const
    {
      return getParallelPieceName( name, path, -1, commSize );
    }

    //! return name of a serial piece file
    /**
     * This is similar to getParallelPieceName, but skips the prefixes for
     * commSize ("s####-") and commRank ("p####-").
     *
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting piece name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     */
    std::string getSerialPieceName(const std::string& name,
                                   const std::string& path) const
    {
      static const std::string extension =
        GridView::dimension == 1 ? ".vtp" : ".vtu";

      return concatPaths(path, name+extension);
    }

    /** \brief write output (interface might change later)
     *
     *  This method can be used in parallel as well as in serial programs.
     *  For serial runs (commSize=1) it chooses other names without the
     *  "s####-p####-" prefix for the .vtu/.vtp files and omits writing of the
     *  .pvtu/pvtp file however.  For parallel runs (commSize > 1) it is the
     *  same as a call to pwrite() with name and path constructed from
     *  a given filename possibly containing a path, and extendpath="".
     *
     *  \param name     Base name of the output files.  This should not
     *                  contain any directory part and no filename extensions.
     *  \param type     How to encode the data in the file.
     *  \param commRank Rank of the current process.
     *  \param commSize Number of processes taking part in this write
     *                  operation.
     */
    std::string write ( const std::string &name,
                        VTK::OutputType type,
                        const int commRank,
                        const int commSize )
    {
      // in the parallel case, just use pwrite, it has all the necessary
      // stuff, so we don't need to reimplement it here.
      if(commSize > 1)
      {
        std::string filename = name;
        std::string path = std::string("");

        // check if a path was already added to name
        // and if yes find filename without path
        auto pos = name.rfind('/');
        if( pos != std::string::npos )
        {
          // extract filename without path
          filename = name.substr( pos+1 );

          // extract the path and added it before
          // the magic below is added
          path = name.substr(0, pos);
        }

        return pwrite(filename, path, "", type, commRank, commSize);
      }

      // make data mode visible to private functions
      outputtype = type;

      // generate filename for process data
      std::string pieceName = getSerialPieceName(name, "");

      // write process data
      std::ofstream file;
      file.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                      std::ios_base::eofbit);
      // check if file can be opened
      try {
        file.open( pieceName.c_str(), std::ios::binary );
      }
      catch(...) {
        std::cerr << "Filename: " << pieceName << " could not be opened" << std::endl;
        throw;
      }
      if (! file.is_open())
        DUNE_THROW(IOError, "Could not write to piece file " << pieceName);
      writeDataFile( file );
      file.close();

      return pieceName;
    }

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
     * \param ot         How to encode the data in the file.
     * \param commRank   Rank of the current process.
     * \param commSize   Number of processes taking part in this write
     *                   operation.
     *
     * \note Currently, extendpath may not be absolute unless path is
     *       absolute, because that would require the value of the current
     *       directory.
     *
     * \throw NotImplemented Extendpath is absolute but path is relative.
     * \throw IOError        Failed to open a file.
     */
    std::string pwrite(const std::string& name, const std::string& path,
                       const std::string& extendpath,
                       VTK::OutputType ot, const int commRank,
                       const int commSize )
    {
      // make data mode visible to private functions
      outputtype=ot;

      // do some magic because paraview can only cope with relative paths to piece files
      std::ofstream file;
      file.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                      std::ios_base::eofbit);
      std::string piecepath = concatPaths(path, extendpath);
      std::string relpiecepath = relativePath(path, piecepath);

      // write this processes .vtu/.vtp piece file
      std::string fullname = getParallelPieceName(name, piecepath, commRank,
                                                  commSize);
      // check if file can be opened
      try {
        file.open(fullname.c_str(),std::ios::binary);
      }
      catch(...) {
        std::cerr << "Filename: " << fullname << " could not be opened" << std::endl;
        throw;
      }
      if (! file.is_open())
        DUNE_THROW(IOError, "Could not write to piecefile file " << fullname);
      writeDataFile(file);
      file.close();
      gridView_.comm().barrier();

      // if we are rank 0, write .pvtu/.pvtp parallel header
      fullname = getParallelHeaderName(name, path, commSize);
      if( commRank  ==0 )
      {
        file.open(fullname.c_str());
        if (! file.is_open())
          DUNE_THROW(IOError, "Could not write to parallel file " << fullname);
        writeParallelHeader(file,name,relpiecepath, commSize );
        file.close();
      }
      gridView_.comm().barrier();
      return fullname;
    }

  private:
    //! write header file in parallel case to stream
    /**
     * Writes a .pvtu/.pvtp file for a collection of concurrently written
     * .vtu/.vtp files.
     *
     * \param s         Stream to write the file contents to.
     * \param piecename Base name of the pieces.  Should not contain a
     *                  directory part or filename extension.
     * \param piecepath Directory part of the pieces.  Since paraview does not
     *                  support absolute paths in parallel collection files,
     *                  this should be a path relative to the directory the
     *                  collection file resides in.  A trailing '/' is
     *                  optional, and an empty value "" is equivalent to the
     *                  value "." except it will look nicer in the collection
     *                  file.
     * \param commSize  Number of processes which are producing the VTK
     *                  output.
     */
    void writeParallelHeader(std::ostream& s, const std::string& piecename,
                             const std::string& piecepath, const int commSize)
    {
      VTK::FileType fileType =
        (n == 1) ? VTK::polyData : VTK::unstructuredGrid;

      VTK::PVTUWriter writer(s, fileType);

      writer.beginMain();

      // PPointData
      {
        std::string scalars, vectors;
        std::tie(scalars,vectors) = getDataNames(vertexdata);
        writer.beginPointData(scalars, vectors);
      }
      for (auto it = vertexdata.begin(),
             end = vertexdata.end();
           it != end;
           ++it)
      {
        unsigned writecomps = it->fieldInfo().size();
        if(writecomps == 2) writecomps = 3;
        writer.addArray(it->name(), writecomps, it->fieldInfo().precision());
      }
      writer.endPointData();

      // PCellData
      {
        std::string scalars, vectors;
        std::tie(scalars,vectors) = getDataNames(celldata);
        writer.beginCellData(scalars, vectors);
      }
      for (auto it = celldata.begin(),
             end = celldata.end();
           it != end;
           ++it)
      {
        unsigned writecomps = it->fieldInfo().size();
        if(writecomps == 2) writecomps = 3;
        writer.addArray(it->name(), writecomps, it->fieldInfo().precision());
      }
      writer.endCellData();

      // PPoints
      writer.beginPoints();
      writer.addArray("Coordinates", 3, coordPrec);
      writer.endPoints();

      // Pieces
      for( int i = 0; i < commSize; ++i )
      {
        const std::string& fullname = getParallelPieceName(piecename,
                                                           piecepath, i,
                                                           commSize);
        writer.addPiece(fullname);
      }

      writer.endMain();
    }

    //! write data file to stream
    void writeDataFile (std::ostream& s)
    {
      VTK::FileType fileType =
        (n == 1) ? VTK::polyData : VTK::unstructuredGrid;

      VTK::VTUWriter writer(s, outputtype, fileType);

      // Grid characteristics
      vertexmapper = new VertexMapper( gridView_, mcmgVertexLayout() );
      if (datamode == VTK::conforming)
      {
        number.resize(vertexmapper->size());
        for (std::vector<int>::size_type i=0; i<number.size(); i++) number[i] = -1;
      }
      countEntities(nvertices, ncells, ncorners);

      writer.beginMain(ncells, nvertices);
      writeAllData(writer);
      writer.endMain();

      // write appended binary data section
      if(writer.beginAppended())
        writeAllData(writer);
      writer.endAppended();

      delete vertexmapper; number.clear();
    }

    void writeAllData(VTK::VTUWriter& writer) {
      // PointData
      writeVertexData(writer);

      // CellData
      writeCellData(writer);

      // Points
      writeGridPoints(writer);

      // Cells
      writeGridCells(writer);
    }

  protected:
    std::string getFormatString() const
    {
      if (outputtype==VTK::ascii)
        return "ascii";
      if (outputtype==VTK::base64)
        return "binary";
      if (outputtype==VTK::appendedraw)
        return "appended";
      if (outputtype==VTK::appendedbase64)
        return "appended";
      DUNE_THROW(IOError, "VTKWriter: unsupported OutputType" << outputtype);
    }

    std::string getTypeString() const
    {
      if (n==1)
        return "PolyData";
      else
        return "UnstructuredGrid";
    }

    //! count the vertices, cells and corners
    virtual void countEntities(int &nvertices_, int &ncells_, int &ncorners_)
    {
      nvertices_ = 0;
      ncells_ = 0;
      ncorners_ = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        ncells_++;
        // because of the use of vertexmapper->map(), this iteration must be
        // in the order of Dune's numbering.
        const int subEntities = it->subEntities(n);
        for (int i=0; i<subEntities; ++i)
        {
          ncorners_++;
          if (datamode == VTK::conforming)
          {
            int alpha = vertexmapper->subIndex(*it,i,n);
            if (number[alpha]<0)
              number[alpha] = nvertices_++;
          }
          else
          {
            nvertices_++;
          }
        }
      }
    }

    template<typename T>
    std::tuple<std::string,std::string> getDataNames(const T& data) const
    {
      std::string scalars = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::scalar)
          {
            scalars = it->name();
            break;
          }

      std::string vectors = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::vector)
          {
            vectors = it->name();
            break;
          }
      return std::make_tuple(scalars,vectors);
    }

    template<typename Data, typename Iterator>
    void writeData(VTK::VTUWriter& writer, const Data& data, const Iterator begin, const Iterator end, int nentries)
    {
      for (auto it = data.begin(),
             iend = data.end();
           it != iend;
           ++it)
      {
        const auto& f = *it;
        VTK::FieldInfo fieldInfo = f.fieldInfo();
        std::size_t writecomps = fieldInfo.size();
        switch (fieldInfo.type())
          {
          case VTK::FieldInfo::Type::scalar:
            break;
          case VTK::FieldInfo::Type::vector:
            // vtk file format: a vector data always should have 3 comps (with
            // 3rd comp = 0 in 2D case)
            if (writecomps > 3)
              DUNE_THROW(IOError,"Cannot write VTK vectors with more than 3 components (components was " << writecomps << ")");
            writecomps = 3;
            break;
          case VTK::FieldInfo::Type::tensor:
            DUNE_THROW(NotImplemented,"VTK output for tensors not implemented yet");
          }
        std::shared_ptr<VTK::DataArrayWriter> p
          (writer.makeArrayWriter(f.name(), writecomps, nentries, fieldInfo.precision()));
        if(!p->writeIsNoop())
          for (Iterator eit = begin; eit!=end; ++eit)
          {
            const Entity & e = *eit;
            f.bind(e);
            f.write(eit.position(),*p);
            f.unbind();
            // vtk file format: a vector data always should have 3 comps
            // (with 3rd comp = 0 in 2D case)
            for (std::size_t j=fieldInfo.size(); j < writecomps; ++j)
              p->write(0.0);
          }
      }
    }

    //! write cell data
    virtual void writeCellData(VTK::VTUWriter& writer)
    {
      if(celldata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(celldata);

      writer.beginCellData(scalars, vectors);
      writeData(writer,celldata,cellBegin(),cellEnd(),ncells);
      writer.endCellData();
    }

    //! write vertex data
    virtual void writeVertexData(VTK::VTUWriter& writer)
    {
      if(vertexdata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(vertexdata);

      writer.beginPointData(scalars, vectors);
      writeData(writer,vertexdata,vertexBegin(),vertexEnd(),nvertices);
      writer.endPointData();
    }

    //! write the positions of vertices
    virtual void writeGridPoints(VTK::VTUWriter& writer)
    {
      writer.beginPoints();

      std::shared_ptr<VTK::DataArrayWriter> p
        (writer.makeArrayWriter("Coordinates", 3, nvertices, coordPrec));
      if(!p->writeIsNoop()) {
        VertexIterator vEnd = vertexEnd();
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {
          int dimw=w;
          for (int j=0; j<std::min(dimw,3); j++)
            p->write((*vit).geometry().corner(vit.localindex())[j]);
          for (int j=std::min(dimw,3); j<3; j++)
            p->write(0.0);
        }
      }
      // free the VTK::DataArrayWriter before touching the stream
      p.reset();

      writer.endPoints();
    }

    //! write the connectivity array
    virtual void writeGridCells(VTK::VTUWriter& writer)
    {
      writer.beginCells();

      // connectivity
      {
        std::shared_ptr<VTK::DataArrayWriter> p1
          (writer.makeArrayWriter("connectivity", 1, ncorners, VTK::Precision::int32));
        if(!p1->writeIsNoop())
          for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
            p1->write(it.id());
      }

      // offsets
      {
        std::shared_ptr<VTK::DataArrayWriter> p2
          (writer.makeArrayWriter("offsets", 1, ncells, VTK::Precision::int32));
        if(!p2->writeIsNoop()) {
          int offset = 0;
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            offset += it->subEntities(n);
            p2->write(offset);
          }
        }
      }

      // types
      if (n>1)
      {
        {
          std::shared_ptr<VTK::DataArrayWriter> p3
            (writer.makeArrayWriter("types", 1, ncells, VTK::Precision::uint8));

          if(!p3->writeIsNoop())
          {
            for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
            {
              int vtktype = VTK::geometryType(it->type());
              p3->write(vtktype);
            }
          }
        }


        // if polyhedron cells found also cell faces need to be written
        if( polyhedralCellsPresent_ )
        {
          writeCellFaces( writer );
        }
      }

      writer.endCells();
    }

  protected:
    bool checkForPolyhedralCells() const
    {
      // check if polyhedron cells are present
      for( const auto& geomType : gridView_.indexSet().types( 0 ) )
      {
        if( VTK::geometryType( geomType ) == VTK::polyhedron )
        {
          return true;
        }
      }
      return false;
    }

    //! write the connectivity array
    virtual void writeCellFaces(VTK::VTUWriter& writer)
    {
      if( ! faceVertices_ )
      {
        faceVertices_.reset( new std::pair< std::vector<int>, std::vector<int> > () );
        // fill face vertex structure
        fillFaceVertices( cornerBegin(), cornerEnd(), gridView_.indexSet(),
                          faceVertices_->first, faceVertices_->second );
      }

      std::vector< int >& faces = faceVertices_->first;
      std::vector< int >& faceOffsets = faceVertices_->second;
      assert( int(faceOffsets.size()) == ncells );

      {
        std::shared_ptr<VTK::DataArrayWriter> p4
          (writer.makeArrayWriter("faces", 1, faces.size(), VTK::Precision::int32));
        if(!p4->writeIsNoop())
        {
          for( const auto& face : faces )
            p4->write( face );
        }
      }

      {
        std::shared_ptr<VTK::DataArrayWriter> p5
          (writer.makeArrayWriter("faceoffsets", 1, ncells, VTK::Precision::int32));
        if(!p5->writeIsNoop())
        {
          for( const auto& offset : faceOffsets )
            p5->write( offset );

          // clear face vertex structure
          faceVertices_.reset();
        }
      }
    }

    template <class CornerIterator, class IndexSet, class T>
    inline void fillFaceVertices( CornerIterator it,
                           const CornerIterator end,
                           const IndexSet& indexSet,
                           std::vector<T>& faces,
                           std::vector<T>& faceOffsets )
    {
      if( n == 3 && it != end )
      {
        // clear output arrays
        faces.clear();
        faces.reserve( 15 * ncells );
        faceOffsets.clear();
        faceOffsets.reserve( ncells );

        int offset = 0;

        Cell element = *it;
        int elIndex = indexSet.index( element );
        std::vector< T > vertices;
        vertices.reserve( 30 );
        for( ; it != end; ++it )
        {
          const Cell& cell = *it ;
          const int cellIndex = indexSet.index( cell ) ;
          if( elIndex != cellIndex )
          {
            fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );

            vertices.clear();
            element = cell ;
            elIndex = cellIndex ;
          }
          vertices.push_back( it.id() );
        }

        // fill faces for last element
        fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );
      }
    }

    template <class Entity, class IndexSet, class T>
    static void fillFacesForElement( const Entity& element,
                                     const IndexSet& indexSet,
                                     const std::vector<T>& vertices,
                                     T& offset,
                                     std::vector<T>& faces,
                                     std::vector<T>& faceOffsets )
    {
      const int dim = n;

      std::map< T, T > vxMap;

      // get number of local faces
      const int nVertices = element.subEntities( dim );
      for( int vx = 0; vx < nVertices; ++ vx )
      {
        const int vxIdx = indexSet.subIndex( element, vx, dim );
        vxMap[ vxIdx ] = vertices[ vx ];
      }

      // get number of local faces
      const int nFaces = element.subEntities( 1 );
      // store number of faces for current element
      faces.push_back( nFaces );
      ++offset;
      // extract each face as a set of vertex indices
      for( int fce = 0; fce < nFaces; ++ fce )
      {
        // obtain face
        const auto face = element.template subEntity< 1 > ( fce );

        // get all vertex indices from current face
        const int nVxFace = face.subEntities( dim );
        faces.push_back( nVxFace );
        ++offset ;
        for( int i=0; i<nVxFace; ++i )
        {
          const T vxIndex = indexSet.subIndex( face, i, dim );
          assert( vxMap.find( vxIndex ) != vxMap.end() );
          faces.push_back( vxMap[ vxIndex ] );
          ++offset ;
        }
      }

      // store face offset for each element
      faceOffsets.push_back( offset );
    }

  protected:
    // the list of registered functions
    std::list<VTKLocalFunction> celldata;
    std::list<VTKLocalFunction> vertexdata;

    // the grid
    GridView gridView_;

    // temporary grid information
    int ncells;
    int nvertices;
    int ncorners;
  private:
    VertexMapper* vertexmapper;
    // in conforming mode, for each vertex id (as obtained by vertexmapper)
    // hold its number in the iteration order (VertexIterator)
    std::vector<int> number;
    VTK::DataMode datamode;
    VTK::Precision coordPrec;

    // true if polyhedral cells are present in the grid
    const bool polyhedralCellsPresent_;

    // pointer holding face vertex connectivity if needed
    std::shared_ptr< std::pair< std::vector<int>, std::vector<int> > > faceVertices_;

  protected:
    VTK::OutputType outputtype;
  };

}

#endif
