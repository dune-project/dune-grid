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

#include <vector>
#include <list>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/path.hh>
#include <dune/common/shared_ptr.hh>
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

   http://www.geophysik.uni-muenchen.de/intranet/it-service/applications/paraview/vtk-file-formats/
   (not available any more)

   http://www.geophysik.uni-muenchen.de/~moder/Paraview/VTK_File_Formats.php
   (alternative)
 */

namespace Dune
{
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
    // extract types
    typedef typename GridView::Grid Grid;
    typedef typename GridView::ctype DT;
    enum { n = GridView::dimension };
    enum { w = GridView::dimensionworld };

    typedef typename GridView::template Codim< 0 >::Entity Cell;
    typedef typename GridView::template Codim< n >::Entity Vertex;
    typedef Cell Entity;

    typedef typename GridView::IndexSet IndexSet;

    static const PartitionIteratorType VTK_Partition = InteriorBorder_Partition;

    typedef typename GridView::template Codim< 0 >
    ::template Partition< VTK_Partition >::Iterator
    GridCellIterator;
    typedef typename GridView::template Codim< n >
    ::template Partition< VTK_Partition >::Iterator
    GridVertexIterator;

    typedef MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout > VertexMapper;

  public:
    typedef Dune::VTKFunction< GridView > VTKFunction;
    typedef shared_ptr< const VTKFunction > VTKFunctionPtr;

  protected:
    typedef typename std::list<VTKFunctionPtr>::const_iterator FunctionIterator;

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
        return GenericReferenceElements<DT,n>::general((*this)->type()).position(0,0);
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
      public ForwardIteratorFacade<VertexIterator, const Entity, const Entity&, int>
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
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++cornerIndexDune;
        const int numCorners = git->template count< n >();
        if( cornerIndexDune == numCorners )
        {
          offset += numCorners;
          cornerIndexDune = 0;

          ++git;
          while( (git != gend) && (git->partitionType() != InteriorEntity) )
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
          visited[vertexmapper.map(*git,cornerIndexDune,n)] = true;
      }
      void increment ()
      {
        switch (datamode)
        {
        case VTK::conforming :
          while(visited[vertexmapper.map(*git,cornerIndexDune,n)])
          {
            basicIncrement();
            if (git == gend) return;
          }
          visited[vertexmapper.map(*git,cornerIndexDune,n)] = true;
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
      const Entity& dereference() const
      {
        return *git;
      }
      //! index of vertex within the entity, in Dune-numbering
      int localindex () const
      {
        return cornerIndexDune;
      }
      //! position of vertex inside the entity
      const FieldVector<DT,n> & position () const
      {
        return GenericReferenceElements<DT,n>::general(git->type())
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
      public ForwardIteratorFacade<CornerIterator, const Entity, const Entity&, int>
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
        const int numCorners = git->template count< n >();
        if( cornerIndexVTK == numCorners )
        {
          offset += numCorners;
          cornerIndexVTK = 0;

          ++git;
          while( (git != gend) && (git->partitionType() != InteriorEntity) )
            ++git;
        }
      }
      bool equals (const CornerIterator & cit) const
      {
        return git == cit.git
               && cornerIndexVTK == cit.cornerIndexVTK
               && datamode == cit.datamode;
      }
      const Entity& dereference() const
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
            number[vertexmapper.map(*git,VTK::renumber(*git,cornerIndexVTK),
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
     */
    explicit VTKWriter ( const GridView &gridView,
                         VTK::DataMode dm = VTK::conforming )
      : gridView_( gridView ),
        datamode( dm )
    { }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p Dune::shared_ptr to the function to visualize
     */
    void addCellData (const VTKFunctionPtr & p)
    {
      celldata.push_back(p);
    }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p The function to visualize.  The VTKWriter object will take
     *          ownership of the VTKFunction *p and delete it when it's done.
     */
    void addCellData (VTKFunction* p)   // DUNE_DEPRECATED
    {
      celldata.push_back(VTKFunctionPtr(p));
    }

    /**
     * @brief Add a grid function (represented by container) that lives on the cells of
     * the grid to the visualization.
     *
     * The container has to have random access via operator[] (e. g. std::vector). The
     * value of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the index (corresponding
     * with the grid view) of the element.
     * For vector valued data all components for an element are assumed to
     * be consecutive.
     *
     * @param v The container with the values of the grid function for each cell.
     * @param name A name to identify the grid function.
     * @param ncomps Number of components (default is 1).
     */
    template<class V>
    void addCellData (const V& v, const std::string &name, int ncomps = 1)
    {
      typedef P0VTKFunction<GridView, V> Function;
      for (int c=0; c<ncomps; ++c) {
        std::stringstream compName;
        compName << name;
        if (ncomps>1)
          compName << "[" << c << "]";
        VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c);
        celldata.push_back(VTKFunctionPtr(p));
      }
    }

    /**
     * @brief Add a grid function that lives on the vertices of the grid to the visualization.
     * @param p The function to visualize.  The VTKWriter object will take
     *          ownership of the VTKFunction *p and delete it when it's done.
     */
    void addVertexData (VTKFunction* p)   // DUNE_DEPRECATED
    {
      vertexdata.push_back(VTKFunctionPtr(p));
    }

    /**
     * @brief Add a grid function that lives on the vertices of the grid to the visualization.
     * @param p Dune::shared_ptr to the function to visualize
     */
    void addVertexData (const VTKFunctionPtr & p)
    {
      vertexdata.push_back(p);
    }

    /**
     * @brief Add a grid function (represented by container) that lives on the vertices of the
     * grid to the visualization output.
     *
     * The container has to have random access via operator[] (e. g. std::vector). The value
     * of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the index (corresponding
     * to the grid view) of the vertex.
     * For vector valued data all components for a vertex are assumed to
     * be consecutive.
     *
     * @param v The container with the values of the grid function for each cell.
     * @param name A name to identify the grid function.
     * @param ncomps Number of components (default is 1).
     */
    template<class V>
    void addVertexData (const V& v, const std::string &name, int ncomps=1)
    {
      typedef P1VTKFunction<GridView, V> Function;
      for (int c=0; c<ncomps; ++c) {
        std::stringstream compName;
        compName << name;
        if (ncomps>1)
          compName << "[" << c << "]";
        VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c);
        vertexdata.push_back(VTKFunctionPtr(p));
      }
    }

    //! clear list of registered functions
    void clear ()
    {
      celldata.clear();
      vertexdata.clear();
    }

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
     *  same as a call to pwrite() with path="" and extendpath="".
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
     *                   realtive to the current directory.
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
    //! return name of a parallel piece file
    /**
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting piece name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     * \param commRank Rank of the process to generate a piece name for.
     * \param commSize Number of processes writing a parallel vtk output.
     */
    std::string getParallelPieceName(const std::string& name,
                                     const std::string& path,
                                     int commRank, int commSize) const
    {
      std::ostringstream s;
      if(path.size() > 0) {
        s << path;
        if(path[path.size()-1] != '/')
          s << '/';
      }
      s << 's' << std::setw(4) << std::setfill('0') << commSize << '-';
      s << 'p' << std::setw(4) << std::setfill('0') << commRank << '-';
      s << name;
      if(GridView::dimension > 1)
        s << ".vtu";
      else
        s << ".vtp";
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
      std::ostringstream s;
      if(path.size() > 0) {
        s << path;
        if(path[path.size()-1] != '/')
          s << '/';
      }
      s << 's' << std::setw(4) << std::setfill('0') << commSize << '-';
      s << name;
      if(GridView::dimension > 1)
        s << ".pvtu";
      else
        s << ".pvtp";
      return s.str();
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
     *  same as a call to pwrite() with path="" and extendpath="".
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
        return pwrite(name, "", "", type, commRank, commSize);

      // make data mode visible to private functions
      outputtype = type;

      // generate filename for process data
      std::string pieceName = getSerialPieceName(name, "");

      // write process data
      std::ofstream file;
      file.open( pieceName.c_str(), std::ios::binary );
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

      // do some magic because paraview can only cope with relative pathes to piece files
      std::ofstream file;
      std::string piecepath = concatPaths(path, extendpath);
      std::string relpiecepath = relativePath(path, piecepath);

      // write this processes .vtu/.vtp piece file
      std::string fullname = getParallelPieceName(name, piecepath, commRank,
                                                  commSize);
      file.open(fullname.c_str(),std::ios::binary);
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
        std::string scalars;
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end();
             ++it)
          if ((*it)->ncomps()==1)
          {
            scalars = (*it)->name();
            break;
          }
        std::string vectors;
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end();
             ++it)
          if ((*it)->ncomps()>1)
          {
            vectors = (*it)->name();
            break;
          }
        writer.beginPointData(scalars, vectors);
      }
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end();
           ++it)
      {
        unsigned writecomps = (*it)->ncomps();
        if(writecomps == 2) writecomps = 3;
        writer.addArray<float>((*it)->name(), writecomps);
      }
      writer.endPointData();

      // PCellData
      {
        std::string scalars;
        for (FunctionIterator it=celldata.begin(); it!=celldata.end();
             ++it)
          if ((*it)->ncomps()==1)
          {
            scalars = (*it)->name();
            break;
          }
        std::string vectors;
        for (FunctionIterator it=celldata.begin(); it!=celldata.end();
             ++it)
          if ((*it)->ncomps()>1)
          {
            vectors = (*it)->name();
            break;
          }
        writer.beginCellData(scalars, vectors);
      }
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it) {
        unsigned writecomps = (*it)->ncomps();
        if(writecomps == 2) writecomps = 3;
        writer.addArray<float>((*it)->name(), writecomps);
      }
      writer.endCellData();

      // PPoints
      writer.beginPoints();
      writer.addArray<float>("Coordinates", 3);
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
      vertexmapper = new VertexMapper( gridView_ );
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
    virtual void countEntities(int &nvertices, int &ncells, int &ncorners)
    {
      nvertices = 0;
      ncells = 0;
      ncorners = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        ncells++;
        // because of the use of vertexmapper->map(), this iteration must be
        // in the order of Dune's numbering.
        for (int i=0; i<it->template count<n>(); ++i)
        {
          ncorners++;
          if (datamode == VTK::conforming)
          {
            int alpha = vertexmapper->map(*it,i,n);
            if (number[alpha]<0)
              number[alpha] = nvertices++;
          }
          else
          {
            nvertices++;
          }
        }
      }
    }

    //! write cell data
    virtual void writeCellData(VTK::VTUWriter& writer)
    {
      if(celldata.size() == 0)
        return;

      std::string scalars = "";
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          scalars = (*it)->name();
          break;
        }
      std::string vectors = "";
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          vectors = (*it)->name();
          break;
        }

      writer.beginCellData(scalars, vectors);
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
      {
        // vtk file format: a vector data always should have 3 comps (with
        // 3rd comp = 0 in 2D case)
        unsigned writecomps = (*it)->ncomps();
        if(writecomps == 2) writecomps = 3;
        shared_ptr<VTK::DataArrayWriter<float> > p
          (writer.makeArrayWriter<float>((*it)->name(), writecomps,
                                         ncells));
        if(!p->writeIsNoop())
          for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
              p->write((*it)->evaluate(j,*i,i.position()));
            // vtk file format: a vector data always should have 3 comps
            // (with 3rd comp = 0 in 2D case)
            for (unsigned j=(*it)->ncomps(); j < writecomps; ++j)
              p->write(0.0);
          }
      }
      writer.endCellData();
    }

    //! write vertex data
    virtual void writeVertexData(VTK::VTUWriter& writer)
    {
      if(vertexdata.size() == 0)
        return;

      std::string scalars = "";
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          scalars = (*it)->name();
          break;
        }
      std::string vectors = "";
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          vectors = (*it)->name();
          break;
        }

      writer.beginPointData(scalars, vectors);
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
      {
        // vtk file format: a vector data always should have 3 comps (with
        // 3rd comp = 0 in 2D case)
        unsigned writecomps = (*it)->ncomps();
        if(writecomps == 2) writecomps = 3;
        shared_ptr<VTK::DataArrayWriter<float> > p
          (writer.makeArrayWriter<float>((*it)->name(), writecomps,
                                         nvertices));
        if(!p->writeIsNoop())
          for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
              p->write((*it)->evaluate(j,*vit,vit.position()));
            // vtk file format: a vector data always should have 3 comps
            // (with 3rd comp = 0 in 2D case)
            for (unsigned j=(*it)->ncomps(); j < writecomps; ++j)
              p->write(0.0);
          }
      }
      writer.endPointData();
    }

    //! write the positions of vertices
    virtual void writeGridPoints(VTK::VTUWriter& writer)
    {
      writer.beginPoints();

      shared_ptr<VTK::DataArrayWriter<float> > p
        (writer.makeArrayWriter<float>("Coordinates", 3, nvertices));
      if(!p->writeIsNoop()) {
        VertexIterator vEnd = vertexEnd();
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {
          int dimw=w;
          for (int j=0; j<std::min(dimw,3); j++)
            p->write(vit->geometry().corner(vit.localindex())[j]);
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
        shared_ptr<VTK::DataArrayWriter<int> > p1
          (writer.makeArrayWriter<int>("connectivity", 1, ncorners));
        if(!p1->writeIsNoop())
          for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
            p1->write(it.id());
      }

      // offsets
      {
        shared_ptr<VTK::DataArrayWriter<int> > p2
          (writer.makeArrayWriter<int>("offsets", 1, ncells));
        if(!p2->writeIsNoop()) {
          int offset = 0;
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            offset += it->template count<n>();
            p2->write(offset);
          }
        }
      }

      // types
      if (n>1)
      {
        shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
          (writer.makeArrayWriter<unsigned char>("types", 1, ncells));
        if(!p3->writeIsNoop())
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            int vtktype = VTK::geometryType(it->type());
            p3->write(vtktype);
          }
      }

      writer.endCells();
    }

  protected:
    // the list of registered functions
    std::list<VTKFunctionPtr> celldata;
    std::list<VTKFunctionPtr> vertexdata;

  private:
    // the grid
    GridView gridView_;

    // temporary grid information
  protected:
    int ncells;
    int nvertices;
    int ncorners;
  private:
    VertexMapper* vertexmapper;
    // in conforming mode, for each vertex id (as obtained by vertexmapper)
    // hold its number in the iteration order (VertexIterator)
    std::vector<int> number;
    VTK::DataMode datamode;
  protected:
    VTK::OutputType outputtype;
  };

}

#endif
