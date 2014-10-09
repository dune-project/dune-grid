// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_VTKSEQUENCE_HH
#define DUNE_VTKSEQUENCE_HH

#include <dune/grid/io/file/vtk/vtksequencewriterbase.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune {

  /**
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://www.vtk.org/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class VTKSequenceWriter :
    public VTKWriter<GridView>,
    public VTKSequenceWriterBase<VTKSequenceWriter<GridView> >
  {
    typedef VTKWriter<GridView> BaseType;
    typedef VTKSequenceWriter<GridView> ThisType;
    friend class VTKSequenceWriterBase<VTKSequenceWriter<GridView> >;
  public:
    using VTKSequenceWriterBase<VTKSequenceWriter<GridView> >::write;
    using BaseType::gridView_;
    using BaseType::getParallelPieceName;
    explicit VTKSequenceWriter ( const GridView &gridView,
                                 const std::string& name,
                                 const std::string& path,
                                 const std::string& extendpath,
                                 VTK::DataMode dm = VTK::conforming )
      : BaseType(gridView,dm),
        VTKSequenceWriterBase<VTKSequenceWriter<GridView> >(name,path,extendpath)
    {}
    ~VTKSequenceWriter() {}
  private:
    std::string getParallelHeaderName(const std::string& name,
                                      const std::string& path) const
    {
      return VTKWriter<GridView>::getParallelHeaderName(name,path,this->gridView_.comm().size());
    }
  };

  /**
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://www.vtk.org/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class SubsamplingVTKSequenceWriter :
    public SubsamplingVTKWriter<GridView>,
    public VTKSequenceWriterBase<SubsamplingVTKSequenceWriter<GridView> >
  {
    typedef SubsamplingVTKWriter<GridView> BaseType;
    typedef SubsamplingVTKSequenceWriter<GridView> ThisType;
    friend class VTKSequenceWriterBase<SubsamplingVTKSequenceWriter<GridView> >;
  public:
    using VTKSequenceWriterBase<SubsamplingVTKSequenceWriter<GridView> >::write;
    using BaseType::gridView_;
    using BaseType::getParallelPieceName;
    explicit SubsamplingVTKSequenceWriter ( const GridView &gridView,
                                            unsigned int level_,
                                            const std::string& name,
                                            const std::string& path,
                                            const std::string& extendpath)
      : BaseType(gridView,level_),
        VTKSequenceWriterBase<SubsamplingVTKSequenceWriter<GridView> >(name,path,extendpath)
    {}
    ~SubsamplingVTKSequenceWriter() {}

  };

} // end namespace Dune

#endif
