// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GMSH_TO_ALU_MAIN_HH
#define DUNE_GRID_UTILITY_GMSH_TO_ALU_MAIN_HH

#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

//! Data handle for load_balance()
template<class Grid>
class RedistributeDataHandle :
  public Dune::CommDataHandleIF<RedistributeDataHandle<Grid>, int>
{
  typedef typename Grid::LocalIdSet IdSet;
  typedef typename IdSet::IdType Id;
  const Grid &grid;
  std::map<Id, int> &elementTags;

public:
  RedistributeDataHandle(const Grid &grid_,
                         std::map<Id, int> &elementTags_) :
    grid(grid_), elementTags(elementTags_)
  { };

  bool contains(int dim, int codim) const
  { return codim == 0; }

  bool fixedsize (int dim, int codim) const { return false; }

  template<class Entity>
  std::size_t size (const Entity &e) const {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (std::size_t " <<
               Dune::className(*this) << "::size(const " <<
               Dune::className<Entity>() << "&) const) should never be "
               "called.");
  }

  std::size_t size (const typename Grid::template Codim<0>::Entity &e) const
  { return elementTags.count(grid.localIdSet().id(e)); }

  template<class MessageBuffer, class Entity>
  void gather(MessageBuffer &buff, const Entity &e) const {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (void "
               << Dune::className(*this) << "::gather(" <<
               Dune::className<MessageBuffer>() << "&, const " <<
               Dune::className<Entity>() << "&) const) should never be "
               "called.");
  }

  template<class MessageBuffer>
  void gather(MessageBuffer &buff,
              const typename Grid::template Codim<0>::Entity &e) const
  {
    typename std::map<Id, int>::const_iterator it =
      elementTags.find(grid.localIdSet().id(e));
    if(it != elementTags.end())
      buff.write(it->second);
  }

  template<class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buff, const Entity &e, std::size_t n) {
    DUNE_THROW(Dune::RangeError, "Nothing needs to be communicated for codim "
               << Entity::codimension << ", so this method (void " <<
               Dune::className(*this) << "::scatter(" <<
               Dune::className<MessageBuffer>() << "&, const " <<
               Dune::className<Entity>() << "&, std::size_t)) should never be "
               "called.");
  }

  template<class MessageBuffer>
  void scatter(MessageBuffer &buff,
               const typename Grid::template Codim<0>::Entity &e,
               std::size_t n)
  {
    switch(n) {
    case 0 : break;
    case 1 : buff.read(elementTags[grid.localIdSet().id(e)]); break;
    default : DUNE_THROW(Dune::RangeError, "At most one data item may be "
                         "communicated!");
    }
  }

  void compress() { }
};

int main(int argc, char **argv) {

  try {
    Dune::MPIHelper::instance(argc, argv);

    //////////////////////////////////////////////////////////////////////
    //
    //  parse commandline
    //

    std::string gmshName;
    std::string aluName;
    std::string tagsName;

    if(argc > 1) gmshName = argv[1];
    if(argc > 2) aluName = argv[2];
    if(argc > 3) tagsName = argv[3];
    if(gmshName == "" || aluName == "" || argc > 4) {
      std::cerr <<
      "Convert a .msh-file into a partitioned ALU-macrogridfile\n"
      "\n"
      "SYNOPSIS:\n"
      "  " << programName << " GMSH_FILENAME ALU_PREFIX [TAGS_PREFIX]\n"
      "\n"
      "PARAMETERS:\n"
      "  GMSH_FILENAME .msh-file to read the grid from.\n"
      "  ALU_PREFIX Filename-prefix for ALUGrid to write the macrogrid to.  Each\n"
      "    process appends its rank to this filename to form something like\n"
      "    ALUPREFIX.RANK.\n"
      "  TAGS_PREFIX If given, read physical entitity numbers from the .msh file\n"
      "    and write them into a series of files of the form TAG_PREFIX.RANK.  The\n"
      "    file format is ASCII with one line per mesh element, each line\n"
      "    consisting of a physical.  The lines are written in level 0 iteration\n"
      "    order for all partition types.\n"
      << std::flush;
      return 1;
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  create Grid
    //

    Dune::GridFactory<Grid> gf;
    std::vector<int> insertionElementTags;

    if(Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
    {
      if(tagsName == "")
        Dune::GmshReader<Grid>::read(gf, gmshName);
      else {
        // make them local since we don't do anything with them anyway atm
        std::vector<int> insertionBoundaryTags;
        Dune::GmshReader<Grid>::read(gf, gmshName, insertionBoundaryTags,
                                     insertionElementTags);
      }
    }

    Dune::shared_ptr<Grid> gridp(gf.createGrid());

    //////////////////////////////////////////////////////////////////////
    //
    //  transfer data into maps
    //

    typedef Grid::LocalIdSet::IdType Id;
    typedef std::map<Id, int> TagMap;
    TagMap elementTagMap;
    if(tagsName != "") {
      typedef Grid::LevelGridView GV;
      const GV &gv = gridp->levelView(0);
      const Grid::LocalIdSet &lis = gridp->localIdSet();
      const GV::Codim<0>::Iterator end = gv.end<0>();
      for(GV::Codim<0>::Iterator it = gv.begin<0>(); it != end; ++it)
        elementTagMap[lis.id(*it)] =
          insertionElementTags[gf.insertionIndex(*it)];
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  actual loadbalancing
    //

    if(tagsName == "")
      gridp->loadBalance();
    else {
      typedef RedistributeDataHandle<Grid> DataHandle;
      DataHandle dh(*gridp, elementTagMap);
      gridp->loadBalance
        (static_cast<Dune::CommDataHandleIF<DataHandle, int>&>(dh));
      // also communicate data on non-interior elements
      gridp->levelView(0).communicate(dh, Dune::InteriorBorder_All_Interface,
                                      Dune::ForwardCommunication);
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  write grid
    //

    {
      std::string path;
      std::string basename;

      std::size_t pos = aluName.rfind('/');
      switch(pos) {
      case std::string::npos :
        path = ".";
        basename = aluName;
        break;
      case 0 :
        path = "/";
        basename = aluName.substr(1);
        break;
      default :
        path = aluName.substr(0, pos);
        basename = aluName.substr(pos+1);
        break;
      }
      gridp->writeMacroGrid(path, basename);
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  write tags
    //

    if(tagsName != "") {
      std::ostringstream s;
      s << tagsName << "." << gridp->comm().rank();
      std::ofstream file(s.str());
      if(!file)
        DUNE_THROW(Dune::IOError, "Can't open tags-file " << s.str() << " for "
                   "writing");

      typedef Grid::LevelGridView GV;
      const GV &gv = gridp->levelView(0);
      const Grid::LocalIdSet &lis = gridp->localIdSet();
      const GV::Codim<0>::Iterator &end = gv.end<0>();
      for(GV::Codim<0>::Iterator it = gv.begin<0>(); it != end; ++it)
        file << elementTagMap[lis.id(*it)] << "\n";

      if(!file)
        DUNE_THROW(Dune::IOError, "Write error while writing tags to "
                   "tags-file " << s.str());
    }
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

#endif // DUNE_GRID_UTILITY_GMSH_TO_ALU_MAIN_HH
