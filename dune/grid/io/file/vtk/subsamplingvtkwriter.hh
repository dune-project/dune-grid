// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
    constexpr static int dim = GridView::dimension;
    constexpr static int dimw = GridView::dimensionworld;
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
     * @param intervals_             A wrapper for the number of refined intervals on one
     *                         axis as returned by either refinementIntervals() or
     *                         refinementLevels().
     * @param coerceToSimplex_ Set this to true to always triangulate elements
     *                         into simplices, even where it's not necessary
     *                         (i.e. for hypercubes).
     * @param coordPrecision   the precision in which to write out coordinates
     *
     * The datamode is always nonconforming.
     */
    explicit SubsamplingVTKWriter (const GridView &gridView,
                                   Dune::RefinementIntervals intervals_, bool coerceToSimplex_ = false,
                                   VTK::Precision coordPrecision = VTK::Precision::float32)
        : Base(gridView, VTK::nonconforming, coordPrecision)
        , intervals(intervals_), coerceToSimplex(coerceToSimplex_)
    {
      if(intervals_.intervals() < 1) {
        DUNE_THROW(Dune::IOError,"SubsamplingVTKWriter: Refinement intervals must be larger than zero! (One interval means no subsampling)");
      }
    }

  private:
    GeometryType subsampledGeometryType(GeometryType geometryType)
    {
      return (geometryType.isCube() && !coerceToSimplex ? geometryType : GeometryTypes::simplex(dim));
    }

    template<typename SubIterator>
    struct IteratorSelector
    {};

    SubElementIterator refinementBegin(const Refinement& refinement, Dune::RefinementIntervals intervals, IteratorSelector<SubElementIterator>)
    {
      return refinement.eBegin(intervals);
    }

    SubVertexIterator refinementBegin(const Refinement& refinement, Dune::RefinementIntervals intervals, IteratorSelector<SubVertexIterator>)
    {
      return refinement.vBegin(intervals);
    }

    SubElementIterator refinementEnd(const Refinement& refinement, Dune::RefinementIntervals intervals, IteratorSelector<SubElementIterator>)
    {
      return refinement.eEnd(intervals);
    }

    SubVertexIterator refinementEnd(const Refinement& refinement, Dune::RefinementIntervals intervals, IteratorSelector<SubVertexIterator>)
    {
      return refinement.vEnd(intervals);
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
        std::shared_ptr<VTK::DataArrayWriter> p
          (writer.makeArrayWriter(f.name(), writecomps, nentries, fieldInfo.precision()));
        if(!p->writeIsNoop())
          for (Iterator eit = begin; eit!=end; ++eit)
          {
            const Entity & e = *eit;
            f.bind(e);
            Refinement &refinement =
              buildRefinement<dim, ctype>(eit->type(),
                                          subsampledGeometryType(eit->type()));
            for(SubIterator sit = refinementBegin(refinement,intervals,sis),
                  send = refinementEnd(refinement,intervals,sis);
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
    virtual void countEntities(int &nvertices_, int &ncells_, int &ncorners_);

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

    Dune::RefinementIntervals intervals;
    bool coerceToSimplex;
  };

  //! count the vertices, cells and corners
  template <class GridView>
  void SubsamplingVTKWriter<GridView>::countEntities(int &nvertices_, int &ncells_, int &ncorners_)
  {
    nvertices_ = 0;
    ncells_ = 0;
    ncorners_ = 0;
    for (CellIterator it=this->cellBegin(); it!=cellEnd(); ++it)
    {
      Refinement &refinement = buildRefinement<dim, ctype>(it->type(), subsampledGeometryType(it->type()));

      ncells_ += refinement.nElements(intervals);
      nvertices_ += refinement.nVertices(intervals);
      ncorners_ += refinement.nElements(intervals) * refinement.eBegin(intervals).vertexIndices().size();
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

    std::shared_ptr<VTK::DataArrayWriter> p
      (writer.makeArrayWriter("Coordinates", 3, nvertices, this->coordPrecision()));
    if(!p->writeIsNoop())
      for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
      {
        Refinement &refinement =
          buildRefinement<dim, ctype>(i->type(),
                                      subsampledGeometryType(i->type()));
        for(SubVertexIterator sit = refinement.vBegin(intervals),
            send = refinement.vEnd(intervals);
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
      std::shared_ptr<VTK::DataArrayWriter> p1
        (writer.makeArrayWriter("connectivity", 1, ncorners, VTK::Precision::int32));
      // The offset within the index numbering
      if(!p1->writeIsNoop()) {
        int offset = 0;
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          GeometryType coercedToType = subsampledGeometryType(i->type());
          Refinement &refinement =
            buildRefinement<dim, ctype>(i->type(), coercedToType);
          for(SubElementIterator sit = refinement.eBegin(intervals),
              send = refinement.eEnd(intervals);
              sit != send; ++sit)
          {
            IndexVector indices = sit.vertexIndices();
            for(unsigned int ii = 0; ii < indices.size(); ++ii)
              p1->write(offset+indices[VTK::renumber(coercedToType, ii)]);
          }
          offset += refinement.nVertices(intervals);
        }
      }
    }

    // offsets
    {
      std::shared_ptr<VTK::DataArrayWriter> p2
        (writer.makeArrayWriter("offsets", 1, ncells,  VTK::Precision::int32));
      if(!p2->writeIsNoop()) {
        // The offset into the connectivity array
        int offset = 0;
        for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
        {
          Refinement &refinement =
            buildRefinement<dim, ctype>(i->type(),
                                        subsampledGeometryType(i->type()));
          unsigned int verticesPerCell =
            refinement.eBegin(intervals).vertexIndices().size();
          for(int element = 0; element < refinement.nElements(intervals);
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
      std::shared_ptr<VTK::DataArrayWriter> p3
        (writer.makeArrayWriter("types", 1, ncells, VTK::Precision::uint8));
      if(!p3->writeIsNoop())
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          GeometryType coerceTo = subsampledGeometryType(it->type());
          Refinement &refinement =
            buildRefinement<dim, ctype>(it->type(), coerceTo);
          int vtktype = VTK::geometryType(coerceTo);
          for(int i = 0; i < refinement.nElements(intervals); ++i)
            p3->write(vtktype);
        }
    }

    writer.endCells();
  }
}

#endif // DUNE_SUBSAMPLINGVTKWRITER_HH
