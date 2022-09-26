// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_VTKSEQUENCEWRITERBASE_HH
#define DUNE_GRID_IO_FILE_VTK_VTKSEQUENCEWRITERBASE_HH

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <memory>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/common/path.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune {

  /** \brief Base class to write pvd-files which contains a list of all collected vtk-files
   *
   * Derive from this class to write pvd-file suitable for easy visualization with
   * <a href="http://www.vtk.org/">The Visualization Toolkit (VTK)</a>.
   *
   * \tparam GridView Grid view of the grid we are writing
   *
   */

  template<class GridView>
  class VTKSequenceWriterBase
  {
    std::shared_ptr<VTKWriter<GridView> > vtkWriter_;
    std::vector<double> timesteps_;
    std::string name_,path_,extendpath_;
    int rank_;
    int size_;
  public:
    /** \brief Set up the VTKSequenceWriterBase class
     *
     * \param vtkWriter Writer object used to write the individual time step data files
     * \param rank Process number in a multi-process setting
     * \param size Total number of processes
     */
    explicit VTKSequenceWriterBase( std::shared_ptr<VTKWriter<GridView> > vtkWriter,
                                    const std::string& name,
                                    const std::string& path,
                                    const std::string& extendpath,
                                    int rank,
                                    int size)
      : vtkWriter_(vtkWriter),
        name_(name), path_(path),
        extendpath_(extendpath),
        rank_(rank),
        size_(size)
    {}

    /**
     * accessor for the underlying VTKWriter instance
     */
    const std::shared_ptr< VTKWriter<GridView> >& vtkWriter() const
    {
      return vtkWriter_;
    }

    /** \brief Adds a field of cell data to the VTK file */
    void addCellData (const std::shared_ptr<const typename VTKWriter<GridView>::VTKFunction> &p)
    {
      vtkWriter_->addCellData(p);
    }

    /** \brief Adds a field of cell data to the VTK file
     * \param v The container with the values of the grid function for each cell
     * \param name A name to identify the grid function
     * \param ncomps Number of components (default is 1)
     */
    template<class V >
    void addCellData (const V &v, const std::string &name, int ncomps=1)
    {
      vtkWriter_->addCellData(v, name, ncomps);
    }

    /** \brief Adds a field of vertex data to the VTK file */
    void addVertexData (const std::shared_ptr<const typename VTKWriter<GridView>::VTKFunction> &p)
    {
      vtkWriter_->addVertexData(p);
    }

    /** \brief Adds a field of vertex data to the VTK file
     * \param v The container with the values of the grid function for each vertex
     * \param name A name to identify the grid function
     * \param ncomps Number of components (default is 1)
     */
    template<class V >
    void addVertexData (const V &v, const std::string &name, int ncomps=1)
    {
      vtkWriter_->addVertexData(v, name, ncomps);
    }


    /**
     * \brief Writes VTK data for the given time,
     * \param time The time(step) for the data to be written.
     * \param type VTK output type.
     */
    void write (double time, VTK::OutputType type = VTK::ascii)
    {
      /* remember current time step */
      unsigned int count = timesteps_.size();
      timesteps_.push_back(time);

      /* write VTK file */
      if(size_==1)
        vtkWriter_->write(concatPaths(path_,seqName(count)),type);
      else
        vtkWriter_->pwrite(seqName(count), path_,extendpath_,type);

      /* write pvd file ... only on rank 0 */
      if (rank_==0) {
        std::ofstream pvdFile;
        pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                           std::ios_base::eofbit);
        std::string pvdname = name_ + ".pvd";
        pvdFile.open(pvdname.c_str());
        pvdFile << "<?xml version=\"1.0\"?> \n"
                << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << VTK::getEndiannessString() << "\"> \n"
                << "<Collection> \n";
        for (unsigned int i=0; i<=count; i++)
        {
          // filename
          std::string piecepath;
          std::string fullname;
          if(size_==1) {
            piecepath = path_;
            fullname = vtkWriter_->getSerialPieceName(seqName(i), piecepath);
          }
          else {
            piecepath = concatPaths(path_, extendpath_);
            fullname = vtkWriter_->getParallelHeaderName(seqName(i), piecepath, size_);
          }
          pvdFile << "<DataSet timestep=\"" << timesteps_[i]
                  << "\" group=\"\" part=\"0\" name=\"\" file=\""
                  << fullname << "\"/> \n";
        }
        pvdFile << "</Collection> \n"
                << "</VTKFile> \n" << std::flush;
        pvdFile.close();
      }
    }

    /**
     * \brief Clears all VTK data added to the VTK writer
     */
    void clear()
    {
      vtkWriter_->clear();
    }

    /**
     * \brief Retrieve the current list of time steps
     */
    const std::vector<double>& getTimeSteps() const
    {
      return timesteps_;
    }

    /**
     * \brief Set the current list of time steps
     * \note This makes it possible to serialize the sequence writers state.
     *       Can be used to continue writing a VTK sequence after a restart of the program.
     */
    void setTimeSteps(const std::vector<double>& timesteps)
    {
      timesteps_ = timesteps;
    }

  private:

    // create sequence name
    std::string seqName(unsigned int count) const
    {
      std::stringstream n;
      n.fill('0');
      n << name_ << "-" << std::setw(5) << count;
      return n.str();
    }
  };

} // end namespace Dune

#endif
