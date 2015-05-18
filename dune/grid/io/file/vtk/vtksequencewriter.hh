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
    public VTKSequenceWriterBase<GridView>
  {
  public:
    explicit VTKSequenceWriter ( const GridView &gridView,
                                 const std::string& name,
                                 const std::string& path,
                                 const std::string& extendpath,
                                 VTK::DataMode dm = VTK::conforming )
      : VTKSequenceWriterBase<GridView>(std::make_shared<VTKWriter<GridView> >(gridView,dm),
                                        name,path,extendpath,
                                        gridView.comm().rank(), gridView.comm().size())
    {}

    ~VTKSequenceWriter() {}
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
    public VTKSequenceWriterBase<GridView>
  {
  public:
    explicit SubsamplingVTKSequenceWriter ( const GridView &gridView,
                                            unsigned int level_,
                                            const std::string& name,
                                            const std::string& path,
                                            const std::string& extendpath)
      : VTKSequenceWriterBase<GridView>(std::make_shared<SubsamplingVTKWriter<GridView> >(gridView,level_),
                                        name,path,extendpath,
                                        gridView.comm().rank(), gridView.comm().size())
    {}
    ~SubsamplingVTKSequenceWriter() {}

  };

} // end namespace Dune

#endif
