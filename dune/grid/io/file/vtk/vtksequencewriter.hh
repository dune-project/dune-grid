// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_VTKSEQUENCE_HH
#define DUNE_VTKSEQUENCE_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>

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
  class VTKSequenceWriter : public VTKWriter<GridView> {
    typedef VTKWriter<GridView> BaseType;
    typedef VTKSequenceWriter<GridView> ThisType;
    std::string name_,path_,extendpath_;
    std::vector<double> timesteps_;
  public:
    explicit VTKSequenceWriter ( const GridView &gridView,
                                 const std::string& name,
                                 const std::string& path,
                                 const std::string& extendpath,
                                 VTK::DataMode dm = VTK::conforming )
      : BaseType(gridView,dm),
        name_(name), path_(path),
        extendpath_(extendpath)
    {}
    ~VTKSequenceWriter() {}

    /**
     * \brief Writes VTK data for the given time.
     * \param time The time(step) for the data to be written.
     * \param ot VTK output type.
     */
    void write (double time, VTK::OutputType ot = VTK::ascii)
    {
      /* remember current time step */
      unsigned int count = timesteps_.size();
      timesteps_.push_back(time);

      /* make sure the directory exists */
      // mkdir("vtk", 777);

      /* write VTK file */
      std::string pvtuName = BaseType::pwrite(seqName(count), path_,extendpath_,ot);

      /* write pvd file ... only on rank 0 */
      if (this->gridView_.comm().rank()==0) {
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
            this->getParallelPieceName(seqName(i), piecepath,
                                       this->gridView_.comm().rank(),
                                       this->gridView_.comm().size());
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

    // do not inherit pwrite
    void pwrite();
  };

} // end namespace Dune

#endif
