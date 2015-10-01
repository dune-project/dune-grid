// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_SUBSAMPLINGVTKWRITER_HH
#define DUNE_SUBSAMPLINGVTKWRITER_HH

#include <ostream>
#include <memory>

#include <dune/common/indent.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/virtualrefinement.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

/** @file
    @author Jö Fahlke
    @brief Provides subsampled file i/o for the visualization toolkit
 */


namespace Dune
{
  /**
   * @brief Writer for the output of subsampled grid functions in the vtk format.
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
    typedef typename GridView::template Codim< 0 >::Entity Entity;
    typedef VirtualRefinement<dim, ctype> Refinement;
    typedef typename Refinement::IndexVector IndexVector;
    typedef typename Refinement::ElementIterator SubElementIterator;
    typedef typename Refinement::VertexIterator SubVertexIterator;

    typedef typename Base::CellIterator CellIterator;
    typedef typename Base::FunctionIterator FunctionIterator;
    using Base::cellBegin;
    using Base::cellEnd;
    using Base::celldata;
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
                                   int level_, bool coerceToSimplex_ = false)
      : Base(gridView, VTK::nonconforming)
        , level(level_), coerceToSimplex(coerceToSimplex_)
    {
      if(level_ < 0) {
        DUNE_THROW(Dune::IOError,"SubsamplingVTKWriter: Negative Subsampling " << level_ << " must not be used!");
      }
    }

  private:
    GeometryType subsampledGeometryType(GeometryType geometryType) {
      if(geometryType.isCube() && !coerceToSimplex) { /* nothing */ }
      else geometryType.makeSimplex(dim);
      return geometryType;
    }


    template<typename SubIterator>
    struct IteratorSelector
    {};

    SubElementIterator refinementBegin(const Refinement& refinement, int level, IteratorSelector<SubElementIterator>)
    {
      return refinement.eBegin(level);
    }

    SubVertexIterator refinementBegin(const Refinement& refinement, int level, IteratorSelector<SubVertexIterator>)
    {
      return refinement.vBegin(level);
    }

    SubElementIterator refinementEnd(const Refinement& refinement, int level, IteratorSelector<SubElementIterator>)
    {
      return refinement.eEnd(level);
    }

    SubVertexIterator refinementEnd(const Refinement& refinement, int level, IteratorSelector<SubVertexIterator>)
    {
      return refinement.vEnd(level);
    }

    template<typename Data, typename Iterator, typename SubIterator>
    void writeData(VTK::VTUWriter& writer, const Data& data, const Iterator begin, const Iterator end, int nentries, IteratorSelector<SubIterator> sis)
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
        std::shared_ptr<VTK::DataArrayWriter<float> > p
          (writer.makeArrayWriter<float>(f.name(), writecomps, nentries));
        if(!p->writeIsNoop())
          for (Iterator eit = begin; eit!=end; ++eit)
          {
            const Entity & e = *eit;
            f.bind(e);
            Refinement &refinement =
              buildRefinement<dim, ctype>(eit->type(),
                                          subsampledGeometryType(eit->type()));
            for(SubIterator sit = refinementBegin(refinement,level,sis),
                  send = refinementEnd(refinement,level,sis);
                sit != send;
                ++sit)
              {
                f.write(sit.coords(),*p);
                // expand 2D-Vectors to 3D for VTK format
                for(unsigned j = f.fieldInfo().size(); j < writecomps; j++)
                  p->write(0.0);
              }
            f.unbind();
          }
      }
    }


  protected:
    //! count the vertices, cells and corners
    virtual void countEntities(int &nvertices, int &ncells, int &ncorners);

    //! write cell data
    virtual void writeCellData(VTK::VTUWriter& writer);

    //! write vertex data
    virtual void writeVertexData(VTK::VTUWriter& writer);

    //! write the positions of vertices
    virtual void writeGridPoints(VTK::VTUWriter& writer);

    //! write the connectivity array
    virtual void writeGridCells(VTK::VTUWriter& writer);

  public:
    using Base::addVertexData;
    using Base::addCellData;

  private:
    // hide addVertexData -- adding raw data directly without a VTKFunction
    // currently does not make sense for subsampled meshes, as the higher order
    // information is missing. See FS#676.
    template<class V>
    void addVertexData (const V& v, const std::string &name, int ncomps=1);
    template<class V>
    void addCellData (const V& v, const std::string &name, int ncomps=1);

    int level;
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
  void SubsamplingVTKWriter<GridView>::writeCellData(VTK::VTUWriter& writer)
  {
    if(celldata.size() == 0)
      return;

    // Find the names of the first scalar and vector data fields.
    // These will be marked as the default fields (the ones that ParaView shows
    // when the file has just been opened).
    std::string defaultScalarField, defaultVectorField;
    std::tie(defaultScalarField, defaultVectorField) = this->getDataNames(celldata);

    writer.beginCellData(defaultScalarField, defaultVectorField);
    writeData(writer,celldata,cellBegin(),cellEnd(),ncells,IteratorSelector<SubElementIterator>());
    writer.endCellData();
  }

  //! write vertex data
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::writeVertexData(VTK::VTUWriter& writer)
  {
    if(vertexdata.size() == 0)
      return;

    // Find the names of the first scalar and vector data fields.
    // These will be marked as the default fields (the ones that ParaView shows
    // when the file has just been opened).
    std::string defaultScalarField, defaultVectorField;
    std::tie(defaultScalarField, defaultVectorField) = this->getDataNames(vertexdata);

    writer.beginPointData(defaultScalarField, defaultVectorField);
    writeData(writer,vertexdata,cellBegin(),cellEnd(),nvertices,IteratorSelector<SubVertexIterator>());
    writer.endPointData();
  }

  //! write the positions of vertices
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::writeGridPoints(VTK::VTUWriter& writer)
  {
    writer.beginPoints();

    std::shared_ptr<VTK::DataArrayWriter<float> > p
      (writer.makeArrayWriter<float>("Coordinates", 3, nvertices));
    if(!p->writeIsNoop())
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement =
          buildRefinement<dim, ctype>(i->type(),
                                      subsampledGeometryType(i->type()));
        for(SubVertexIterator sit = refinement.vBegin(level),
            send = refinement.vEnd(level);
            sit != send; ++sit)
        {
          FieldVector<ctype, dimw> coords = i->geometry().global(sit.coords());
          for (int j=0; j<std::min(int(dimw),3); j++)
            p->write(coords[j]);
          for (int j=std::min(int(dimw),3); j<3; j++)
            p->write(0.0);
        }
      }
    // free the VTK::DataArrayWriter before touching the stream
    p.reset();

    writer.endPoints();
  }

  //! write the connectivity array
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::writeGridCells(VTK::VTUWriter& writer)
  {
    writer.beginCells();

    // connectivity
    {
      std::shared_ptr<VTK::DataArrayWriter<int> > p1
        (writer.makeArrayWriter<int>("connectivity", 1, ncorners));
      // The offset within the index numbering
      if(!p1->writeIsNoop()) {
        int offset = 0;
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          GeometryType coercedToType = subsampledGeometryType(i->type());
          Refinement &refinement =
            buildRefinement<dim, ctype>(i->type(), coercedToType);
          for(SubElementIterator sit = refinement.eBegin(level),
              send = refinement.eEnd(level);
              sit != send; ++sit)
          {
            IndexVector indices = sit.vertexIndices();
            for(unsigned int ii = 0; ii < indices.size(); ++ii)
              p1->write(offset+indices[VTK::renumber(coercedToType, ii)]);
          }
          offset += refinement.nVertices(level);
        }
      }
    }

    // offsets
    {
      std::shared_ptr<VTK::DataArrayWriter<int> > p2
        (writer.makeArrayWriter<int>("offsets", 1, ncells));
      if(!p2->writeIsNoop()) {
        // The offset into the connectivity array
        int offset = 0;
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          Refinement &refinement =
            buildRefinement<dim, ctype>(i->type(),
                                        subsampledGeometryType(i->type()));
          unsigned int verticesPerCell =
            refinement.eBegin(level).vertexIndices().size();
          for(int element = 0; element < refinement.nElements(level);
              ++element)
          {
            offset += verticesPerCell;
            p2->write(offset);
          }
        }
      }
    }

    // types
    if (dim>1)
    {
      std::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
        (writer.makeArrayWriter<unsigned char>("types", 1, ncells));
      if(!p3->writeIsNoop())
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          GeometryType coerceTo = subsampledGeometryType(it->type());
          Refinement &refinement =
            buildRefinement<dim, ctype>(it->type(), coerceTo);
          int vtktype = VTK::geometryType(coerceTo);
          for(int i = 0; i < refinement.nElements(level); ++i)
            p3->write(vtktype);
        }
    }

    writer.endCells();
  }
}

#endif // DUNE_SUBSAMPLINGVTKWRITER_HH
