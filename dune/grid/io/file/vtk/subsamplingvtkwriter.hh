// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_SUBSAMPLINGVTKWRITER_HH
#define DUNE_SUBSAMPLINGVTKWRITER_HH

#include <ostream>

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/virtualrefinement.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

/** @file
    @author Jö Fahlke
    @brief Provides subsampled file i/o for the visualization toolkit
 */


namespace Dune
{
  /**
   * @brief Writer for the ouput of subsampled grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit
   * (VTK)</a>.  In contrast to the regular VTKWriter, this Writer allows
   * subsampling of the elements via VirtualRefinement.  The
   * SubSamplingVTKWriter always writes nonconforming data.
   */
  template< class GridView >
  class SubsamplingVTKWriter
    : public VTKWriter<GridView>
  {
    typedef VTKWriter<GridView> Base;
    enum { dim = GridView::dimension };
    enum { dimw = GridView::dimensionworld };
    typedef typename GridView::Grid::ctype ctype;
    typedef VirtualRefinement<dim, ctype> Refinement;
    typedef typename Refinement::IndexVector IndexVector;
    typedef typename Refinement::ElementIterator SubElementIterator;
    typedef typename Refinement::VertexIterator SubVertexIterator;

    typedef typename Base::CellIterator CellIterator;
    typedef typename Base::FunctionIterator FunctionIterator;
    typedef typename Base::SimpleStream SimpleStream;
    using Base::cellBegin;
    using Base::cellEnd;
    using Base::celldata;
    using Base::indent;
    using Base::indentDown;
    using Base::indentUp;
    using Base::ncells;
    using Base::ncorners;
    using Base::nvertices;
    using Base::outputtype;
    using Base::vertexBegin;
    using Base::vertexEnd;
    using Base::vertexdata;

  public:
    /**
     * @brief Construct a SubsamplingVTKWriter working on a specific GridView.
     *
     * @param gridView         The gridView the grid functions live
     *                         on. (E. g. a LevelGridView.)
     * @param level_           The level for the subrefinement.
     * @param coerceToSimplex_ Set this to true to always triangulate elements
     *                         into simplices, even where it's not necessary
     *                         (i.e. for hypercubes).
     *
     * The datamode is always nonconforming.
     */
    explicit SubsamplingVTKWriter (const GridView &gridView,
                                   unsigned int level_, bool coerceToSimplex_ = false)
      : Base(gridView, VTK::nonconforming)
        , level(level_), coerceToSimplex(coerceToSimplex_)
    { }

  private:
    GeometryType subsampledGeometryType(GeometryType geometryType) {
      if(geometryType.isCube() && !coerceToSimplex) { /* nothing */ }
      else geometryType.makeSimplex(dim);
      return geometryType;
    }

  protected:
    //! count the vertices, cells and corners
    virtual void countEntities(int &nvertices, int &ncells, int &ncorners);

    //! write cell data
    virtual void writeCellData(std::ostream& s,
                               VTKDataArrayWriterFactory& factory);

    //! write vertex data
    virtual void writeVertexData(std::ostream& s,
                                 VTKDataArrayWriterFactory& factory);

    //! write the positions of vertices
    virtual void writeGridPoints(std::ostream& s,
                                 VTKDataArrayWriterFactory& factory);

    //! write the connectivity array
    virtual void writeGridCells(std::ostream& s,
                                VTKDataArrayWriterFactory& factory);

    //! write the appended data sections
    virtual void writeAppendedData (std::ostream& s);

  public:
    using Base::addVertexData;

  private:
    // hide addVertexData -- adding vertex data directly without a VTKFunction
    // currently does not work since the P1VectorWrapper used for that uses a
    // nearest-neighbour search to find the value for the given point.  See
    // FS#676.
    template<class V>
    void addVertexData (const V& v, const std::string &name, int ncomps=1);

    unsigned int level;
    bool coerceToSimplex;
  };

  //! count the vertices, cells and corners
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::countEntities(int &nvertices, int &ncells, int &ncorners)
  {
    nvertices = 0;
    ncells = 0;
    ncorners = 0;
    for (CellIterator it=this->cellBegin(); it!=cellEnd(); ++it)
    {
      Refinement &refinement = buildRefinement<dim, ctype>(it->type(), subsampledGeometryType(it->type()));

      ncells += refinement.nElements(level);
      nvertices += refinement.nVertices(level);
      ncorners += refinement.nElements(level) * refinement.eBegin(level).vertexIndices().size();
    }
  }

  //! write cell data
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::
  writeCellData(std::ostream& s, VTKDataArrayWriterFactory& factory)
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
    s << ">" << std::endl;
    indentUp();
    for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
    {
      // vtk file format: a vector data always should have 3 comps (with 3rd
      // comp = 0 in 2D case)
      unsigned writecomps = (*it)->ncomps();
      if(writecomps == 2) writecomps = 3;

      shared_ptr<VTKDataArrayWriter<float> > p
        (factory.make<float>((*it)->name(), writecomps, ncells));
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        for(SubElementIterator sit = refinement.eBegin(level), send = refinement.eEnd(level);
            sit != send; ++sit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            p->write((*it)->evaluate(j,*i,sit.coords()));
          // expand 2D-Vectors to 3D
          for(unsigned j = (*it)->ncomps(); j < writecomps; j++)
            p->write(0.0);
        }
      }
    }
    indentDown();
    indent(s); s << "</CellData>" << std::endl;
  }

  //! write vertex data
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::
  writeVertexData(std::ostream& s, VTKDataArrayWriterFactory& factory)
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
    s << ">" << std::endl;
    indentUp();
    for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
    {
      // vtk file format: a vector data always should have 3 comps (with 3rd
      // comp = 0 in 2D case)
      unsigned writecomps = (*it)->ncomps();
      if(writecomps == 2) writecomps = 3;

      shared_ptr<VTKDataArrayWriter<float> > p
        (factory.make<float>((*it)->name(), writecomps, nvertices));
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        for(SubVertexIterator sit = refinement.vBegin(level), send = refinement.vEnd(level);
            sit != send; ++sit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            p->write((*it)->evaluate(j,*i,sit.coords()));
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          for(unsigned j = (*it)->ncomps(); j < writecomps; j++)
            p->write(0.0);
        }
      }
    }
    indentDown();
    indent(s); s << "</PointData>" << std::endl;
  }

  //! write the positions of vertices
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::
  writeGridPoints(std::ostream& s, VTKDataArrayWriterFactory& factory)
  {
    indent(s); s << "<Points>" << std::endl;
    indentUp();

    shared_ptr<VTKDataArrayWriter<float> > p
      (factory.make<float>("Coordinates", 3, nvertices));
    for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
    {
      Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
      for(SubVertexIterator sit = refinement.vBegin(level), send = refinement.vEnd(level);
          sit != send; ++sit)
      {
        FieldVector<ctype, dimw> coords = i->geometry().global(sit.coords());
        for (int j=0; j<std::min(int(dimw),3); j++)
          p->write(coords[j]);
        for (int j=std::min(int(dimw),3); j<3; j++)
          p->write(0.0);
      }
    }
    // free the VTKDataArrayWriter before touching the stream
    p.reset();

    indentDown();
    indent(s); s << "</Points>" << std::endl;
  }

  //! write the connectivity array
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::
  writeGridCells(std::ostream& s, VTKDataArrayWriterFactory& factory)
  {
    indent(s);
    if (dim>1)
      s << "<Cells>" << std::endl;
    else
      s << "<Lines>" << std::endl;
    indentUp();

    // connectivity
    {
      shared_ptr<VTKDataArrayWriter<int> > p1
        (factory.make<int>("connectivity", 1, ncorners));
      // The offset within the index numbering
      int offset = 0;
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        GeometryType coercedToType = subsampledGeometryType(i->type());
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), coercedToType);
        for(SubElementIterator sit = refinement.eBegin(level), send = refinement.eEnd(level);
            sit != send; ++sit)
        {
          IndexVector indices = sit.vertexIndices();
          for(unsigned int ii = 0; ii < indices.size(); ++ii)
            p1->write(offset+indices[VTK::renumber(coercedToType, ii)]);
        }
        offset += refinement.nVertices(level);
      }
    }

    // offsets
    {
      shared_ptr<VTKDataArrayWriter<int> > p2
        (factory.make<int>("offsets", 1, ncells));
      // The offset into the connectivity array
      int offset = 0;
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        unsigned int verticesPerCell = refinement.eBegin(level).vertexIndices().size();
        for(int element = 0; element < refinement.nElements(level); ++element)
        {
          offset += verticesPerCell;
          p2->write(offset);
        }
      }
    }

    // types
    if (dim>1)
    {
      shared_ptr<VTKDataArrayWriter<unsigned char> > p3
        (factory.make<unsigned char>("types", 1, ncells));
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        GeometryType coerceTo = subsampledGeometryType(it->type());
        Refinement &refinement = buildRefinement<dim, ctype>(it->type(), coerceTo);
        int vtktype = VTK::geometryType(coerceTo);
        for(int i = 0; i < refinement.nElements(level); ++i)
          p3->write(vtktype);
      }
    }

    indentDown();
    indent(s);
    if (dim>1)
      s << "</Cells>" << std::endl;
    else
      s << "</Lines>" << std::endl;
  }

  //! write the appended data sections
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::writeAppendedData (std::ostream& s)
  {
    indent(s); s << "<AppendedData encoding=\"raw\">" << std::endl;
    indentUp();
    indent(s); s << "_"; // indicates start of binary data

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
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        for(SubVertexIterator sit = refinement.vBegin(level), send = refinement.vEnd(level);
            sit != send; ++sit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            stream.write(float((*it)->evaluate(j,*i,sit.coords())));
          //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if((*it)->ncomps()==2)
            stream.write(float(0.0));
        }
      }
    }

    // cell data
    for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
    {
      blocklength = ncells * (*it)->ncomps() * sizeof(float);
      //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
      if((*it)->ncomps()==2)
        blocklength = ncells * (3) * sizeof(float);
      stream.write(blocklength);
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        for(SubElementIterator sit = refinement.eBegin(level), send = refinement.eEnd(level);
            sit != send; ++sit)
        {
          for (int j=0; j<(*it)->ncomps(); j++)
            stream.write(float((*it)->evaluate(j,*i,sit.coords())));
          if((*it)->ncomps() == 2)
            stream.write(float(0.0));
        }
      }
    }

    // point coordinates
    blocklength = nvertices * 3 * sizeof(float);
    stream.write(blocklength);
    for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
    {
      Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
      for(SubVertexIterator sit = refinement.vBegin(level), send = refinement.vEnd(level);
          sit != send; ++sit)
      {
        FieldVector<ctype, dimw> coords = i->geometry().global(sit.coords());
        for (int j=0; j<std::min(int(dimw),3); j++)
          stream.write(float(coords[j]));
        for (int j=std::min(int(dimw),3); j<3; j++)
          stream.write(float(0.0));
      }
    }

    // connectivity
    blocklength = ncorners * sizeof(unsigned int);
    stream.write(blocklength);
    {
      //offset into the vertex numbering
      unsigned int offset = 0;
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        GeometryType coerceTo = subsampledGeometryType(i->type());
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), coerceTo);
        for(SubElementIterator sit = refinement.eBegin(level), send = refinement.eEnd(level);
            sit != send; ++sit)
        {
          IndexVector indices = sit.vertexIndices();
          for(unsigned int ii = 0; ii < indices.size(); ++ii)
            stream.write(offset+indices[VTK::renumber(coerceTo, ii)]);
        }
        offset += refinement.nVertices(level);
      }
    }

    // offsets
    blocklength = ncells * sizeof(unsigned int);
    stream.write(blocklength);
    {
      //offset into the connectivity array
      unsigned int offset = 0;
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement = buildRefinement<dim, ctype>(i->type(), subsampledGeometryType(i->type()));
        unsigned int verticesPerCell = refinement.eBegin(level).vertexIndices().size();
        for(int element = 0; element < refinement.nElements(level); ++element)
        {
          offset += verticesPerCell;
          stream.write(offset);
        }
      }
    }

    // cell types
    if (dim>1)
    {
      blocklength = ncells * sizeof(unsigned char);
      stream.write(blocklength);
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        GeometryType coerceTo = subsampledGeometryType(it->type());
        Refinement &refinement = buildRefinement<dim, ctype>(it->type(), coerceTo);
        unsigned char vtktype = VTK::geometryType(coerceTo);
        for(int i = 0; i < refinement.nElements(level); ++i)
          stream.write(vtktype);
      }
    }

    s << std::endl;
    indentDown();
    indent(s); s << "</AppendedData>" << std::endl;
  }
}

#endif // DUNE_SUBSAMPLINGVTKWRITER_HH
