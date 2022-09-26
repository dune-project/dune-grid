// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_VTKSEQUENCE_HH
#define DUNE_VTKSEQUENCE_HH

#include <memory>

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
    /** \brief Constructor with a given VTKWriter or SubsamplingVTKWriter
     *
     * At each time step, the VTKSequenceWriter writes the grid and data currently attached to the vtkWriter object.
     * All calls to the addCellData and addVertexData methods of the VTKSequenceWriter class are forwarded to the
     * vtkWriter, but we propose that you call the corresponding methods on the vtkWriter directly.
     *
     * \param name       Base name of the output files.  This should not
     *                   contain any directory part and not filename
     *                   extensions.  It will be used both for each processes
     *                   piece as well as the parallel collection file.
     */
    VTKSequenceWriter ( std::shared_ptr<VTKWriter<GridView> > vtkWriter,
                        const std::string& name )
      : VTKSequenceWriterBase<GridView>(vtkWriter,
                                        name,
                                        "",
                                        "",
                                        vtkWriter->gridView_.comm().rank(),
                                        vtkWriter->gridView_.comm().size())
    {}

    /** \brief Constructor with a given VTKWriter or SubsamplingVTKWriter
     *
     * At each time step, the VTKSequenceWriter writes the grid and data currently attached to the vtkWriter object.
     * All calls to the addCellData and addVertexData methods of the VTKSequenceWriter class are forwarded to the
     * vtkWriter, but we propose that you call the corresponding methods on the vtkWriter directly.
     *
     * \param name       Base name of the output files.  This should not
     *                   contain any directory part and not filename
     *                   extensions.  It will be used both for each processes
     *                   piece as well as the parallel collection file.
     * \param path       Directory where to put the parallel collection
     *                   (.pvtu/.pvtp) files.  If it is relative, it is taken
     *                   relative to the current directory.
     * \param extendpath Directory where to put the piece files (.vtu/.vtp) of
     *                   this process.  If it is relative, it is taken
     *                   relative to the directory denoted by path.
     */
    VTKSequenceWriter ( std::shared_ptr<VTKWriter<GridView> > vtkWriter,
                        const std::string& name,
                        const std::string& path,
                        const std::string& extendpath )
      : VTKSequenceWriterBase<GridView>(vtkWriter,
                                        name,
                                        path,
                                        extendpath,
                                        vtkWriter->gridView_.comm().rank(),
                                        vtkWriter->gridView_.comm().size())
    {}

    /** \brief Constructor creating its own VTKWriter object
     *
     * At each time step, the VTKSequenceWriter writes the grid and data currently attached to the vtkWriter object.
     * All calls to the addCellData and addVertexData methods of the VTKSequenceWriter class are forwarded to the
     * vtkWriter.
     */
    explicit VTKSequenceWriter ( const GridView &gridView,
                                 const std::string& name,
                                 const std::string& path,
                                 const std::string& extendpath,
                                 VTK::DataMode dm = VTK::conforming )
      : VTKSequenceWriterBase<GridView>(std::make_shared<VTKWriter<GridView> >(gridView,dm),
                                        name,path,extendpath,
                                        gridView.comm().rank(), gridView.comm().size())
    {}
  };

} // end namespace Dune

#endif
