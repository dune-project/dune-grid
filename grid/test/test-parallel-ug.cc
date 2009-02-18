// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: test-ug.cc 4424 2008-09-29 07:46:41Z sander $
// Test parallel interface if a parallel UG is used

#include <config.h>

#include <iostream>
#include <vector>

#include <dune/grid/uggrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/mpihelper.hh>

using namespace Dune;

//! Use only elements for a mapper
template <int dim>
struct ElementLayout
{
  bool contains (Dune::GeometryType gt)
  {
    return gt.dim() == dim;
  }
};

template <int dim>
struct NodeLayout
{
  bool contains (Dune::GeometryType gt)
  {
    return gt.dim() == 0;
  }
};

// A DataHandle class to exchange entries of a vector.
template<class MapperT, int commCodim>
class DataExchange
  : public Dune::CommDataHandleIF<DataExchange<MapperT,
            commCodim>,
        double>
{
public:
  //! export type of data for message buffer
  typedef double DataType;
  typedef std::vector<Dune::FieldVector<double,1> > UserDataType;

  //! constructor
  DataExchange(const MapperT &mapper,
               UserDataType &userDataSend,
               UserDataType &userDataReceive)
    : mapper_(mapper),
      userDataSend_(userDataSend),
      userDataReceive_(userDataReceive)
  {}


  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return (codim == commCodim);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    return true;
  }

  /*! how many objects of type DataType have to be sent for a given entity

     Note: Only the sender side needs to know this size.
   */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    return 1;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather(MessageBuffer& buff, const EntityType& e) const
  {
    DataType x = mapper_.map(e);

    userDataSend_[mapper_.map(e)][0] = x;
    std::cout << "Process "
              << Dune::MPIHelper::getCollectiveCommunication().rank()+1
              << " sends for entity "
              << mapper_.map(e)
              << ": "
              << std::setprecision(20)
              << x << "\n";

    buff.write(x);
  }

  /*! unpack data from message buffer to user

     n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class EntityType>
  void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
  {
    DataType x;
    buff.read(x);

    userDataReceive_[mapper_.map(e)][0] = x;
    std::cout << "Process "
              << Dune::MPIHelper::getCollectiveCommunication().rank()+1
              << " received for entity "
              << mapper_.map(e)
              << ": "
              << std::setprecision(20)
              << x << "\n";
  }

private:
  const MapperT &mapper_;
  UserDataType &userDataSend_;
  UserDataType &userDataReceive_;
};

int main (int argc , char **argv) try
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper &mpiHelper =
    Dune::MPIHelper::instance(argc, argv);

  int rank = mpiHelper.rank();
  int size = mpiHelper.size();

  std::cout << "This is process " << rank << " of " << size << ", PID " << getpid() << " .\n";

  ////////////////////////////////////////////////////////////
  // Create uniform grid on rank 0
  ////////////////////////////////////////////////////////////
  int n = 11;
  const int dim = 2;
  const int commCodim = 0;    // codim of the entities to communicate
  typedef UGGrid<dim> GridType;
  GridType grid;
  GridFactory<GridType> factory(&grid);

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      FieldVector<double,dim> pos;
      pos[0] = double(i)/(n-1);
      pos[1] = double(j)/(n-1);
      factory.insertVertex(pos);
    }
  }

  for (int i=0; i<n-1; i++) {
    for (int j=0; j<n-1; j++) {
      std::vector<unsigned int> v(4);
      v[0] = i*n + j;
      v[1] = i*n + j+1;
      v[2] = (i+1)*n + j;
      v[3] = (i+1)*n + j+1;
      factory.insertElement(GeometryType(GeometryType::cube,dim), v);
    }
  }

  factory.createGrid();


  //////////////////////////////////////////////////////
  // Distribute the grid
  //////////////////////////////////////////////////////
  /*    grid.loadBalance(0,   // strategy
                       0,   // minlevel
                       2,   // depth
                       32,   // maxlevel
                       1);   // minelement
   */
  grid.loadBalance();


  std::cout << "Process " << rank << " has " << grid.size(0) << " elements." << std::endl;

  // make sure each process has some elements
  assert(grid.size(0) > 0);
  // make sure each process has some nodes
  assert(grid.size(dim) > 0);

  //////////////////////////////////////////////////////
  // Test level and leaf mappers/gridViews
  //////////////////////////////////////////////////////
  GridType::LevelGridView level0GridView = grid.levelView(0);
  GridType::LeafGridView leafGridView = grid.leafView();

  std::cout << "LevelGridView for level 0 has " << level0GridView.size(0)
            << " elements and " << level0GridView.size(dim) << " nodes.\n";
  std::cout << "LeafGridView has " << leafGridView.size(0)
            << " elements and " << leafGridView.size(dim) << " nodes.\n";

  typedef LevelMultipleCodimMultipleGeomTypeMapper<GridType,ElementLayout > LevelMapperType;
  LevelMapperType levelMapper(grid, 0);

  typedef LeafMultipleCodimMultipleGeomTypeMapper<GridType,ElementLayout > LeafMapperType;
  LeafMapperType leafMapper(grid);

  std::cout << "Level index set has " << levelMapper.size() << " entities\n";
  std::cout << "Leaf index set has " << leafMapper.size() << " entities\n";

  // create the user data arrays
  typedef std::vector<Dune::FieldVector<double, 1> > UserDataType;
  UserDataType userDataSend(leafGridView.size(commCodim), 0.0);
  UserDataType userDataReceive(leafGridView.size(commCodim), 0.0);
  UserDataType partitionType(leafGridView.size(commCodim), -1e10);

  // write the partition type of each node into the corrosponding
  // result array
  GridType::LeafGridView::Codim<commCodim>::Iterator
    it = leafGridView.begin<commCodim>();
  const GridType::LeafGridView::Codim<commCodim>::Iterator
  &endIt = leafGridView.end<commCodim>();
  for (; it != endIt; ++it) {
    partitionType[leafMapper.map(*it)] = it->partitionType();
  };

  //////////////////////////////////////////////////////
  // Test communication
  //////////////////////////////////////////////////////

  // initialize data handle (marks the nodes where some data was
  // send or received)
  typedef DataExchange<LeafMapperType, commCodim> MyDataHandle;
  MyDataHandle datahandle(leafMapper,
                          userDataSend,
                          userDataReceive);
  // communicate the entities at the interior border to all other
  // processes
  grid.communicate(datahandle,
                   Dune::InteriorBorder_All_Interface,
                   Dune::ForwardCommunication);

  //////////////////////////////////////////////////////
  // Write results to disc
  //////////////////////////////////////////////////////
  std::cout << "Writing data to disk\n";
  Dune::VTKWriter<GridType::LeafGridView> writer(leafGridView);
  if (commCodim == 0) {
    writer.addCellData(userDataSend, "Send");
    writer.addCellData(userDataReceive, "Receive");
    writer.addCellData(partitionType, "Partition Type");
  }
  else if (commCodim == dim) {
    writer.addVertexData(userDataSend, "Send");
    writer.addVertexData(userDataReceive, "Receive");
    writer.addVertexData(partitionType, "Partition Type");
  }
  writer.write("test-parallel-ug",
               Dune::VTKOptions::ascii);
  std::cout << "Done writing data to disk\n";

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
