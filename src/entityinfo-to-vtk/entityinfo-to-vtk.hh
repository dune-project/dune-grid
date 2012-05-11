// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_HH
#define DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_HH

#include <cstddef>
#include <string>
#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

inline int partitionTypeToInt(Dune::PartitionType pt) {
  switch(pt) {
  case Dune::InteriorEntity : return 0;
  case Dune::BorderEntity :   return 1;
  case Dune::OverlapEntity :  return 2;
  case Dune::FrontEntity :    return 3;
  case Dune::GhostEntity :    return 4;
  default :                   return -1;
  }
}

template<Dune::PartitionIteratorType Partition, class GV>
void entityinfoToVTK(const GV &gv, const std::string &vtkName,
                     Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw)
{

  //////////////////////////////////////////////////////////////////////
  //
  //  collect information
  //

  typedef Dune::MultipleCodimMultipleGeomTypeMapper
  <GV, Dune::MCMGElementLayout> EMapper;
  EMapper eMapper(gv);
  std::vector<int> eRanks(eMapper.size(), gv.comm().rank());
  std::vector<int> ePartitionTypes(eMapper.size());
  std::vector<typename GV::IndexSet::IndexType> eIndices(eMapper.size());
  std::vector<int> eMap(eMapper.size());
  std::vector<unsigned> eTopologyIds(eMapper.size());
  {
    typedef typename GV::template Codim<0>::
    template Partition<Partition>::Iterator Iterator;
    const Iterator &end = gv.template end<0, Partition>();
    for(Iterator it = gv.template begin<0, Partition>(); it != end; ++it) {
      int map = eMapper.map(*it);
      ePartitionTypes[map] = partitionTypeToInt(it->partitionType());
      eIndices[map] = gv.indexSet().index(*it);
      eMap[map] = map;
      eTopologyIds[map] = it->type().id();
    }
  }

  typedef Dune::MultipleCodimMultipleGeomTypeMapper
  <GV, Dune::MCMGVertexLayout> VMapper;
  VMapper vMapper(gv);
  std::vector<int> vRanks(vMapper.size(), gv.comm().rank());
  std::vector<int> vPartitionTypes(vMapper.size());
  std::vector<typename GV::IndexSet::IndexType> vIndices(vMapper.size());
  std::vector<int> vMap(vMapper.size());
  {
    static const std::size_t dim = GV::dimension;
    typedef typename GV::template Codim<dim>::
    template Partition<Partition>::Iterator Iterator;
    const Iterator &end = gv.template end<dim, Partition>();
    for(Iterator it = gv.template begin<dim, Partition>(); it != end; ++it) {
      int map = vMapper.map(*it);
      vPartitionTypes[map] = partitionTypeToInt(it->partitionType());
      vIndices[map] = gv.indexSet().index(*it);
      vMap[map] = map;
    }
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

  Dune::VTKWriter<GV, Partition> vtkWriter(gv);

  vtkWriter.addCellData(eRanks, "eRanks");
  vtkWriter.addCellData(ePartitionTypes, "ePartitionTypes");
  vtkWriter.addCellData(eIndices, "eIndices");
  vtkWriter.addCellData(eMap, "eMap");
  vtkWriter.addCellData(eTopologyIds, "eTopologyIds");

  vtkWriter.addVertexData(vRanks, "vRanks");
  vtkWriter.addVertexData(vPartitionTypes, "vPartitionTypes");
  vtkWriter.addVertexData(vIndices, "vIndices");
  vtkWriter.addVertexData(vMap, "vMap");

  vtkWriter.pwrite(basename, path, "", outputtype);
}

template<class GV>
void entityinfoToVTK(const GV &gv, const std::string &vtkName,
                     Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw)
{
  entityinfoToVTK<Dune::All_Partition>(gv, vtkName, outputtype);
}

#endif // DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_HH
