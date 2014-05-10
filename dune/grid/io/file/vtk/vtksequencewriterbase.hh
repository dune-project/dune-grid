// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_VTKSEQUENCEBASE_HH
#define DUNE_VTKSEQUENCEBASE_HH

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/common/path.hh>

namespace Dune {

  /** \brief Base class to write pvd-files which contains a list of all collected vtk-files
   *
   * Derive from this class to write pvd-file suitable for easy visualization with
   * <a href="http://www.vtk.org/">The Visualization Toolkit (VTK)</a>.
   * The derived class needs to inherit from the VTKWriter or the SubsamplingVTKWriter
   *
   * \tparam Imp Type of the derived class (CRTP-trick)
   *
   */

  template<class Imp>
  class VTKSequenceWriterBase
  {
    std::vector<double> timesteps_;
    std::string name_,path_,extendpath_;
  public:
    explicit VTKSequenceWriterBase( const std::string& name,
                                    const std::string& path,
                                    const std::string& extendpath)
      : name_(name), path_(path),
        extendpath_(extendpath)
    {}
    ~VTKSequenceWriterBase() {}

    /**
     * \brief Writes VTK data for the given time,
     * \param time The time(step) for the data to be written.
     * \param ot VTK output type.
     */
    void write (double time, VTK::OutputType ot = VTK::ascii)
    {
      /* remember current time step */
      unsigned int count = timesteps_.size();
      timesteps_.push_back(time);

      /* write VTK file */
      std::string pvtuName = static_cast<Imp*>(this)->pwrite(seqName(count), path_,extendpath_,ot);

      /* write pvd file ... only on rank 0 */
      if (static_cast<Imp*>(this)->gridView_.comm().rank()==0) {
        std::ofstream pvdFile;
        pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                           std::ios_base::eofbit);
        std::string pvdname = name_ + ".pvd";
        pvdFile.open(pvdname.c_str());
        pvdFile << "<?xml version=\"1.0\"?> \n"
                << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"> \n"
                << "<Collection> \n";
        for (unsigned int i=0; i<=count; i++)
        {
          // filename
          std::string piecepath = concatPaths(path_, extendpath_);
          std::string fullname =
            static_cast<Imp*>(this)->getParallelPieceName(seqName(i), piecepath,
                                                          static_cast<Imp*>(this)->gridView_.comm().rank(),
                                                          static_cast<Imp*>(this)->gridView_.comm().size());
          pvdFile << "<DataSet timestep=\"" << timesteps_[i]
                  << "\" group=\"\" part=\"0\" name=\"\" file=\""
                  << fullname << "\"/> \n";
        }
        pvdFile << "</Collection> \n"
                << "</VTKFile> \n" << std::flush;
        pvdFile.close();
      }
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
