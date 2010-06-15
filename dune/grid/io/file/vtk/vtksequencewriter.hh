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
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class VTKSequenceWriter : public VTKWriter<GridView> {
    typedef VTKWriter<GridView> BaseType;
    typedef VTKSequenceWriter<GridView> ThisType;
    std::string name_,path_,extendpath_;
    unsigned int count_;
    std::ofstream seqFile_;
  public:
    explicit VTKSequenceWriter ( const GridView &gridView,
                                 const std::string& name,
                                 const std::string& path,
                                 const std::string& extendpath,
                                 VTK::DataMode dm = VTK::conforming )
      : BaseType(gridView,dm),
        name_(name), path_(path),
        extendpath_(extendpath),
        count_(0),
        seqFile_()
    {
      if (gridView.comm().rank()==0) {
        std::stringstream pvdname;
        pvdname << name << ".pvd";
        seqFile_.open(pvdname.str().c_str());
        seqFile_ << "<?xml version=\"1.0\"?> \n"
                 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"> \n"
                 << "<Collection> \n" << std::flush;
      }
    }
    ~VTKSequenceWriter() {
      if (seqFile_.is_open()) {
        seqFile_ << "</Collection> \n"
                 << "</VTKFile> \n" << std::flush;
      }
    }
    void write (double time, VTK::OutputType ot = VTK::ascii)
    {
      std::stringstream name;
      name.fill('0');
      name << name_ << "-" << std::setw(5) << count_;
      std::string pvtuName = BaseType::pwrite(name.str().c_str(),
                                              path_.c_str(),extendpath_.c_str(),ot);
      if (seqFile_.is_open()) {
        seqFile_ << "<DataSet timestep=\"" << time
                 << "\" group=\"\" part=\"0\" name=\"\" file=\""
                 << pvtuName << "\"/> \n" << std::flush;
      }
      ++count_;
    }
  private:
    // do not inherite pwrite
    void pwrite();
  };
}
#endif
