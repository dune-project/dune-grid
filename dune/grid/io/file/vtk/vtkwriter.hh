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
#include <dune/common/iteratorfacades.hh>
#include <dune/common/path.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
#include <dune/grid/io/file/vtk/function.hh>

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
    typedef typename Grid::ctype DT;
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
    typedef Dune::VTKFunction<Grid> VTKFunction;
    typedef shared_ptr< Dune::VTKFunction<Grid> > VTKFunctionPtr;

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
      CellIterator(const GridCellIterator & x) : GridCellIterator(x) {};
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
      public ForwardIteratorFacade<VertexIterator, Entity, Entity&, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTKOptions::DataMode datamode;
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
                     const VTKOptions::DataMode & dm,
                     const VertexMapper & vm) :
        git(x), gend(end), datamode(dm), cornerIndexDune(0),
        vertexmapper(vm), visited(vm.size(), false),
        offset(0)
      {
        if (datamode == VTKOptions::conforming && git != gend)
          visited[vertexmapper.map(*git,cornerIndexDune,n)] = true;
      };
      void increment ()
      {
        switch (datamode)
        {
        case VTKOptions::conforming :
          while(visited[vertexmapper.map(*git,cornerIndexDune,n)])
          {
            basicIncrement();
            if (git == gend) return;
          }
          visited[vertexmapper.map(*git,cornerIndexDune,n)] = true;
          break;
        case VTKOptions::nonconforming :
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
      Entity& dereference() const
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
      public ForwardIteratorFacade<CornerIterator, Entity, Entity&, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTKOptions::DataMode datamode;
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
                     const VTKOptions::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), cornerIndexVTK(0),
        vertexmapper(vm),
        number(num), offset(0) {};
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
      Entity& dereference() const
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
        case VTKOptions::conforming :
          return
            number[vertexmapper.map(*git,renumber(*git,cornerIndexVTK),n)];
        case VTKOptions::nonconforming :
          return offset + renumber(*git,cornerIndexVTK);
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
                         VTKOptions::DataMode dm = VTKOptions::conforming )
      : gridView_( gridView ),
        datamode( dm )
    {
      indentCount = 0;
      numPerLine = 4*3; //should be a multiple of 3 !
    }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p Dune:shared_ptr to the function to visualize
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
     * @param p Dune:shared_ptr to the function to visualize
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
     *  "s####:p####:" prefix for the .vtu/.vtp files and omits writing of the
     *  .pvtu/pvtp file however.  For parallel runs (commSize > 1) it is the
     *  same as a call to pwrite() with path="" and extendpath="".
     *
     *  \param[in]  name  basic name to write (may not contain a path)
     *  \param[in]  type  type of output (e.g,, ASCII) (optional)
     */
    std::string write ( const std::string &name,
                        VTKOptions::OutputType type = VTKOptions::ascii )
    {
      return write( name, type, gridView_.comm().rank(), gridView_.comm().size() );
    }

    /** \brief write output (interface might change later)
     *
     * "pwrite" means "path write" (i.e. write somewhere else than the current
     * directory).  The "p" does not mean this method has a monopoly on
     * parallel writing, the regular write(const std::string &,
     * VTKOptions::OutputType) method ca do that just fine.
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
    std::string pwrite ( const char* name,  const char* path, const char* extendpath,
                         VTKOptions::OutputType type = VTKOptions::ascii )
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
      s << 's' << std::setw(4) << std::setfill('0') << commSize << ':';
      s << 'p' << std::setw(4) << std::setfill('0') << commRank << ':';
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
      s << 's' << std::setw(4) << std::setfill('0') << commSize << ':';
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
     * commSize ("s####:") and commRank ("p####:").
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
     *  "s####:p####:" prefix for the .vtu/.vtp files and omits writing of the
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
                        VTKOptions::OutputType type,
                        const int commRank,
                        const int commSize )
    {
      // in the parallel case, just use pwrite, it has all the necessary
      // stuff, so we don't need to reimplement it here.
      if(commSize > 1)
        return pwrite(name, "", "", type, commRank, commSize);

      // make data mode visible to private functions
      outputtype = type;

      // reset byte counter for binary appended output
      bytecount = 0;

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
                       VTKOptions::OutputType ot, const int commRank,
                       const int commSize )
    {
      // make data mode visible to private functions
      outputtype=ot;

      // reset byte counter for binary appended output
      bytecount = 0;

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
      // xml header
      s << "<?xml version=\"1.0\"?>\n";

      // VTKFile
      s << "<VTKFile type=\"P" << getTypeString()
        << "\" version=\"0.1\" byte_order=\"" << getEndiannessString() << "\">\n";
      indentUp();

      // PUnstructuredGrid
      indent(s);
      s << "<P" << getTypeString() << " GhostLevel=\"0\">\n";
      indentUp();

      // PPointData
      indent(s); s << "<PPointData";
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          s << " Scalars=\"" << (*it)->name() << "\"" ;
          break;
        }
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          s << " Vectors=\"" << (*it)->name() << "\"" ;
          break;
        }
      s << ">\n";
      indentUp();
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
      {
        indent(s); s << "<PDataArray type=\"Float32\" Name=\"" << (*it)->name() << "\" ";
        s << "NumberOfComponents=\"" << ((*it)->ncomps()>1 ? 3 : 1) << "\" ";
        s << "format=\"" << getFormatString() << "\"/>\n";
      }
      indentDown();
      indent(s); s << "</PPointData>\n";

      // PCellData
      indent(s); s << "<PCellData";
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          s << " Scalars=\"" << (*it)->name() << "\"" ;
          break;
        }
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          s << " Vectors=\"" << (*it)->name() << "\"" ;
          break;
        }
      s << ">\n";
      indentUp();
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
      {
        indent(s); s << "<PDataArray type=\"Float32\" Name=\"" << (*it)->name() << "\" ";
        s << "NumberOfComponents=\"" << ((*it)->ncomps()>1 ? 3 : 1) << "\" ";
        s << "format=\"" << getFormatString() << "\"/>\n";
      }
      indentDown();
      indent(s); s << "</PCellData>\n";

      // PPoints
      indent(s); s << "<PPoints>\n";
      indentUp();
      indent(s); s << "<PDataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"" << "3" << "\" ";
      s << "format=\"" << getFormatString() << "\"/>\n";
      indentDown();
      indent(s); s << "</PPoints>\n";

      // Pieces
      for( int i = 0; i < commSize; ++i )
      {
        const std::string& fullname = getParallelPieceName(piecename,
                                                           piecepath, i,
                                                           commSize);
        indent(s); s << "<Piece Source=\"" << fullname << "\"/>\n";
      }

      // /PUnstructuredGrid
      indentDown();
      indent(s);
      s << "</P" << getTypeString() << ">\n";

      // /VTKFile
      indentDown();
      s << "</VTKFile>\n";

      s.flush();
    }

    //! write data file to stream
    void writeDataFile (std::ostream& s)
    {
      // xml header
      s << "<?xml version=\"1.0\"?>\n";

      // VTKFile
      s << "<VTKFile type=\"" << getTypeString()
        << "\" version=\"0.1\" byte_order=\"" << getEndiannessString() << "\">\n";
      indentUp();

      // Grid characteristics
      vertexmapper = new VertexMapper( gridView_ );
      if (datamode == VTKOptions::conforming)
      {
        number.resize(vertexmapper->size());
        for (std::vector<int>::size_type i=0; i<number.size(); i++) number[i] = -1;
      }
      countEntities(nvertices, ncells, ncorners);

      // UnstructuredGrid
      indent(s);
      s << "<" << getTypeString() << ">\n";
      indentUp();

      // Piece
      indent(s);
      if (n>1)
        s << "<Piece NumberOfPoints=\"" << nvertices << "\" NumberOfCells=\"" << ncells << "\">\n";
      else
        s << "<Piece NumberOfPoints=\"" << nvertices << "\""
          << " NumberOfVerts=\"0\""
          << " NumberOfLines=\"" << ncells << "\">"
          << " NumberOfPolys=\"0\"\n";
      indentUp();

      // PointData
      writeVertexData(s);

      // CellData
      writeCellData(s);

      // Points
      writeGridPoints(s);

      // Cells
      writeGridCells(s);

      // /Piece
      indentDown();
      indent(s); s << "</Piece>\n";

      // /UnstructuredGrid
      indentDown();
      indent(s);
      s << "</" << getTypeString() << ">\n";

      // write appended binary dat section
      if (outputtype==VTKOptions::binaryappended)
        writeAppendedData(s);

      // /VTKFile
      indentDown();
      s << "</VTKFile>\n";

      s.flush();

      delete vertexmapper; number.clear();
    }

  protected:
    std::string getEndiannessString() const
    {
      short i = 1;
      if (reinterpret_cast<char*>(&i)[1] == 1)
        return "BigEndian";
      else
        return "LittleEndian";
    }

    std::string getFormatString() const
    {
      if (outputtype==VTKOptions::ascii)
        return "ascii";
      if (outputtype==VTKOptions::binary)
        return "binary";
      if (outputtype==VTKOptions::binaryappended)
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
          if (datamode == VTKOptions::conforming)
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
    virtual void writeCellData (std::ostream& s)
    {
      indent(s); s << "<CellData";
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          s << " Scalars=\"" << (*it)->name() << "\"" ;
          break;
        }
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          s << " Vectors=\"" << (*it)->name() << "\"" ;
          break;
        }
      s << ">\n";
      indentUp();
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
      {
        VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, (*it)->name().c_str(), (*it)->ncomps(), (*it)->ncomps()*ncells);
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            p->write((*it)->evaluate(j,*i,i.position()));
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if((*it)->ncomps()==2)
            p->write(0.0);
        }
        delete p;
      }
      indentDown();
      indent(s); s << "</CellData>\n";
      s.flush();
    }

    //! write vertex data
    virtual void writeVertexData (std::ostream& s)
    {
      indent(s); s << "<PointData";
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()==1)
        {
          s << " Scalars=\"" << (*it)->name() << "\"" ;
          break;
        }
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        if ((*it)->ncomps()>1)
        {
          s << " Vectors=\"" << (*it)->name() << "\"" ;
          break;
        }
      s << ">\n";
      indentUp();
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
      {
        VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, (*it)->name().c_str(), (*it)->ncomps(), (*it)->ncomps()*nvertices);
        for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            p->write((*it)->evaluate(j,*vit,vit.position()));
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if((*it)->ncomps()==2)
            p->write(0.0);
        }
        delete p;
      }
      indentDown();
      indent(s); s << "</PointData>\n";
      s.flush();
    }

    //! write the positions of vertices
    virtual void writeGridPoints (std::ostream& s)
    {
      indent(s); s << "<Points>\n";
      indentUp();

      VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, "Coordinates", 3, 3*nvertices);
      VertexIterator vEnd = vertexEnd();
      for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
      {
        int dimw=w;
        for (int j=0; j<std::min(dimw,3); j++)
          p->write(vit->geometry().corner(vit.localindex())[j]);
        for (int j=std::min(dimw,3); j<3; j++)
          p->write(0.0);
      }
      delete p;

      indentDown();
      indent(s); s << "</Points>\n";
      s.flush();
    }

    //! write the connectivity array
    virtual void writeGridCells (std::ostream& s)
    {
      indent(s);
      if (n>1)
        s << "<Cells>\n";
      else
        s << "<Lines>\n";
      indentUp();

      // connectivity
      VTKDataArrayWriter<int> *p1 = makeVTKDataArrayWriter<int>(s, "connectivity", 1, ncorners);
      for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
        p1->write(it.id());
      delete p1;

      // offsets
      VTKDataArrayWriter<int> *p2 = makeVTKDataArrayWriter<int>(s, "offsets", 1, ncells);
      {
        int offset = 0;
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          offset += it->template count<n>();
          p2->write(offset);
        }
      }
      delete p2;

      // types
      if (n>1)
      {
        VTKDataArrayWriter<unsigned char> *p3 = makeVTKDataArrayWriter<unsigned char>(s, "types", 1, ncells);
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          int vtktype = vtkType(it->type());
          p3->write(vtktype);
        }
        delete p3;
      }

      indentDown();
      indent(s);
      if (n>1)
        s << "</Cells>\n";
      else
        s << "</Lines>\n";
      s.flush();
    }

    //! write the appended data sections
    virtual void writeAppendedData (std::ostream& s)
    {
      indent(s); s << "<AppendedData encoding=\"raw\">\n";
      indentUp();
      indent(s); s << "_";   // indicates start of binary data

      SimpleStream stream(s);

      // write length before each data block
      unsigned int blocklength;

      // point data
      for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
      {

        blocklength = nvertices * (*it)->ncomps() * sizeof(float);
        //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
        if((*it)->ncomps()==2)
          blocklength = nvertices * (3) * sizeof(float);
        stream.write(blocklength);
        std::vector<bool> visited(vertexmapper->size(), false);
        for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
          {
            float data = (*it)->evaluate(j,*vit,vit.position());
            stream.write(data);
          }
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if((*it)->ncomps()==2)
          {
            float data=0.0;
            stream.write(data);
          }
        }
      }

      // cell data
      for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
      {
        blocklength = ncells * (*it)->ncomps() * sizeof(float);
        // vtk file format: a vector data always should have 3 comps (with
        // 3rd comp = 0 in 2D case)
        if((*it)->ncomps()==2)
          blocklength = ncells * (3) * sizeof(float);
        stream.write(blocklength);
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
          {
            float data = (*it)->evaluate(j,*i,i.position());
            stream.write(data);
          }
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if((*it)->ncomps()==2)
          {
            float data=0.0;
            stream.write(data);
          }
        }
      }

      // point coordinates
      blocklength = nvertices * 3 * sizeof(float);
      stream.write(blocklength);
      std::vector<bool> visited(vertexmapper->size(), false);
      for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
      {
        int dimw=w;
        float data;
        for (int j=0; j<std::min(dimw,3); j++)
        {
          data = vit->geometry().corner(vit.localindex())[j];
          stream.write(data);
        }
        data = 0;
        for (int j=std::min(dimw,3); j<3; j++)
          stream.write(data);
      }

      // connectivity
      blocklength = ncorners * sizeof(unsigned int);
      stream.write(blocklength);
      for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
      {
        stream.write(it.id());
      }

      // offsets
      blocklength = ncells * sizeof(unsigned int);
      stream.write(blocklength);
      {
        int offset = 0;
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          offset += it->template count<n>();
          stream.write(offset);
        }
      }

      // cell types
      if (n>1)
      {
        blocklength = ncells * sizeof(unsigned char);
        stream.write(blocklength);
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          unsigned char vtktype = vtkType(it->type());
          stream.write(vtktype);
        }
      }

      s << std::endl;
      indentDown();
      indent(s); s << "</AppendedData>\n";
      s.flush();
    }

  protected:
    /** @brief Make a VTKDataArrayWriter with new
     *
     * @param s           The stream to write to
     * @param name        The name of the vtk array
     * @param components  The number of components of the vector
     * @param totallength the total number of entries, i.e. components*vectors
     */
    template<class T>
    VTKDataArrayWriter<T> *makeVTKDataArrayWriter(std::ostream &s,
                                                  const char *name,
                                                  unsigned int components,
                                                  unsigned int totallength)
    {
      return Dune::makeVTKDataArrayWriter<T>(outputtype, s, name, components,
                                             totallength, bytecount);
    }

    //! write out data in binary
    class SimpleStream
    {
    public:
      //! make a new stream
      SimpleStream (std::ostream& theStream)
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

    //! increase indentation level
    void indentUp ()
    {
      indentCount++;
    }

    //! decrease indentation level
    void indentDown ()
    {
      indentCount--;
    }

    //! write indentation to stream
    void indent (std::ostream& s)
    {
      for (int i=0; i<indentCount; i++)
        s << "  ";
    }

    //! renumber VTK <-> Dune
    /**
     * Since the renumbering never does anything more complex than exchanging
     * two indices, this method works both ways.
     */
    static int renumber (const GeometryType &t, int i)
    {
      static const int quadRenumbering[4] = {0,1,3,2};
      static const int cubeRenumbering[8] = {0,1,3,2,4,5,7,6};
      static const int prismRenumbering[6] = {0,2,1,3,5,4};
      static const int pyramidRenumbering[5] = {0,1,3,2,4};
      if (t.isQuadrilateral())
        return quadRenumbering[i];
      if (t.isPyramid())
        return pyramidRenumbering[i];
      if (t.isPrism())
        return prismRenumbering[i];
      if (t.isHexahedron())
        return cubeRenumbering[i];
      return i;
    }
    static int renumber (const Entity& e, int i)
    {
      return renumber(e.type(), i);
    }

    // the list of registered functions
    std::list<VTKFunctionPtr> celldata;
    std::list<VTKFunctionPtr> vertexdata;

  private:
    // the grid
    GridView gridView_;

    // indend counter
    int indentCount;
    int numPerLine;

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
    VTKOptions::DataMode datamode;
  protected:
    VTKOptions::OutputType outputtype;
  private:
    unsigned int bytecount;
  };

  /** \brief VTKWriter on the leaf grid
      \ingroup VTK
   */
  template< class Grid >
  class LeafVTKWriter
    : public VTKWriter< typename Grid::LeafGridView >
  {
    typedef VTKWriter< typename Grid::LeafGridView > Base;

  public:
    /** \brief Construct a VTK writer for the leaf level of a given grid */
    explicit LeafVTKWriter ( const Grid &grid,
                             VTKOptions::DataMode dm = VTKOptions::conforming ) DUNE_DEPRECATED
      : Base( grid.leafView(), dm )
    {}
  };

  /** \brief VTKWriter on a given level grid
      \ingroup VTK
   */
  template< class Grid >
  class LevelVTKWriter
    : public VTKWriter< typename Grid::LevelGridView >
  {
    typedef VTKWriter< typename Grid::LevelGridView > Base;

  public:
    /** \brief Construct a VTK writer for a certain level of a given grid */
    LevelVTKWriter ( const Grid &grid, int level,
                     VTKOptions::DataMode dm = VTKOptions::conforming ) DUNE_DEPRECATED
      : Base( grid.levelView( level ), dm )
    {}
  };

}

#endif
