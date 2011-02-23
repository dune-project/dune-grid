// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_ALU_TO_VTK_MAIN_VOLUME_HH
#define DUNE_GRID_UTILITY_ALU_TO_VTK_MAIN_VOLUME_HH

#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

int main(int argc, char **argv) {

  try {
    Dune::MPIHelper::instance(argc, argv);

    //////////////////////////////////////////////////////////////////////
    //
    //  parse commandline
    //

    std::string aluName;
    std::string vtkName;
    std::string tagsName;

    if(argc > 1) aluName = argv[1];
    if(argc > 2) vtkName = argv[2];
    if(argc > 3) tagsName = argv[3];
    if(aluName == "" || vtkName == "" || argc > 4) {
      std::cerr <<
      "Convert a ALU-macrogridfile to VTK\n"
      "\n"
      "SYNOPSIS:\n"
      "  " << programName << " ALU_PREFIX VTK_PREFIX [TAGS_PREFIX]\n"
      "\n"
      "PARAMETERS:\n"
      "  ALU_PREFIX Filename-prefix for ALUGrid to read the macrogrid from.  Each\n"
      "    process appends its rank to this filename to form something like\n"
      "    ALUPREFIX.RANK.\n"
      "  VTK_FILENAME Filename-prefix for VTK to write to.\n"
      "  TAGS_PREFIX If given, read element tag numbers from a series of files of\n"
      "    the form TAG_PREFIX.RANK.  The file format is the same the output of\n"
      "    gmsh-to-alu: ASCII with one line per mesh element, each line consisting\n"
      "    of a whitespace-separated pair of local-ID and physical entity number.\n"
      << std::flush;
      return 1;
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  create Grid
    //

    std::ostringstream rankName;
    rankName << aluName << "."
             << Dune::MPIHelper::getCollectiveCommunication().rank();

    Grid grid(rankName.str());
    typedef Grid::LeafGridView GV;
    const GV &gv = grid.leafView();

    //////////////////////////////////////////////////////////////////////
    //
    //  read tags
    //

    std::vector<int> elementTags;
    if(tagsName != "") {
      std::map<Grid::LocalIdSet::IdType, int> tagMap;

      std::ostringstream s;
      s << tagsName << "." << grid.comm().rank();
      std::ifstream file(s.str());
      if(!file)
        DUNE_THROW(Dune::IOError, s.str() << ": Can't open tags-file for "
                   "reading");
      for(std::size_t lineNo = 0; !file.eof(); ++lineNo) {
        std::string buf;
        std::getline(file, buf);
        // Don't complain on conversion errors (fail()) at eof, since we will
        // get a conversion error for the
        // zero-length-line-without-final-newline at eof.
        if(file.bad() || (file.fail() && !file.eof()))
          DUNE_THROW(Dune::IOError, s.str() << ":" << lineNo << ": Read "
                     "error.");

        // skip empty lines and comments
        std::size_t pos = buf.find_first_not_of(" \t");
        if(pos == std::string::npos || buf[pos] == '#')
          continue;

        Grid::LocalIdSet::IdType id;
        int tag;
        std::istringstream line(buf);
        line >> id >> tag;
        bool fail = !line;
        // make sure there is no garbage at eol except whitespace: try to
        // extract another char which should fail.
        if(!fail) {
          char c;
          line >> c;
          fail = line;
        }
        if(fail)
          DUNE_THROW(Dune::IOError, s.str() << ":" << lineNo << ": Invalid "
                     "formatted input line (expected \"localId tag\\n\")");
        tagMap[id] = tag;

      }

      // copy into vector
      const GV::Codim<0>::Iterator &end = gv.end<0>();
      const Grid::LocalIdSet &lid = grid.localIdSet();
      typedef Dune::MultipleCodimMultipleGeomTypeMapper<
          GV, Dune::MCMGElementLayout
          > Mapper;
      Mapper mapper(gv);
      elementTags.resize(mapper.size(), 0);
      for(GV::Codim<0>::Iterator it = gv.begin<0>(); it != end; ++it)
        elementTags[mapper.map(*it)] = tagMap[lid.id(*it)];

      // no need to communicate ghost data -- ghosts won't be written anyway
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  write grid
    //

    std::string path;
    std::string basename;

    std::size_t pos = vtkName.rfind('/');
    switch(pos) {
    case std::string::npos :
      path = "";
      basename = vtkName;
      break;
    case 0 :
      path = "/";
      basename = vtkName.substr(1);
      break;
    default :
      path = vtkName.substr(0, pos);
      basename = vtkName.substr(pos+1);
      break;
    }

    Dune::VTKWriter<GV> vtkWriter(gv);
    if(tagsName != "")
      vtkWriter.addCellData(elementTags, "elementTags");
    vtkWriter.pwrite(basename, path, "", Dune::VTK::appendedraw);
  }
  catch(const Dune::Exception &e) {
    std::cerr << "Caught Dune exception: " << e << std::endl;
    throw;
  }
  catch(const std::exception &e) {
    std::cerr << "Caught std::exception: " << e.what() << std::endl;
    throw;
  }
  catch(...) {
    std::cerr << "Caught unknown exception" << std::endl;
    throw;
  }
}

#endif // DUNE_GRID_UTILITY_ALU_TO_VTK_MAIN_VOLUME_HH
