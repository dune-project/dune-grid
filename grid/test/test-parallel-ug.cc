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

//! a "flexible" layout where for which arbitrary codims can be
//! specified
template <int commCodim>
struct LayoutWrapper
{
  template <int dim>
  struct Layout
  {
    bool contains(Dune::GeometryType gt)
    { return gt.dim() == dim - commCodim;  }
  };
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


template <class GridView>
void checkConsistency(const GridView &gv)
{
  typedef typename GridView::template Codim<0>::Entity Element;
  typedef typename GridView::template Codim<0>::Iterator Iterator;
  typedef typename Element::LeafIntersectionIterator IntersectionIterator;
  Iterator it = gv.template begin<0>();
  const Iterator &endIt = gv.template end<0>();
  for (; it != endIt; ++it) {
    IntersectionIterator isIt = it->ileafbegin();
    const IntersectionIterator &isEndIt = it->ileafend();
    int n = 0;
    for (; isIt != isEndIt; ++isIt) {
      isIt->boundary();
      isIt->inside();
      if (isIt->neighbor()) {
        isIt->outside();
      }
      ++ n;
    };

    assert(n == 2*GridView::dimension);     // quadrilaterals
  }
};

template <class Grid>
void insertVertices(int n, Dune::GridFactory<Grid> &factory, int coordIdx=0)
{
  static FieldVector<double,Grid::dimension> pos;
  if (coordIdx == Grid::dimension) {
    factory.insertVertex(pos);
    return;
  }

  for (int i=0; i < n; ++i) {
    pos[coordIdx] = double(i)/(n-1);
    insertVertices(n, factory, coordIdx + 1);
  }
}

template <class GridType>
void insertElements(int n, Dune::GridFactory<GridType> &factory, unsigned coordIdx=0)
{
  const int dim = GridType::dimension;

  if (dim == 3) {
    for (int i=0; i<n-1; i++) {
      for (int j=0; j<n-1; j++) {
        for (int k=0; k<n-1; k++) {
          std::vector<unsigned int> v(8);
          v[0] = k*n*n + i*n + j;
          v[1] = k*n*n + i*n + j+1;
          v[2] = k*n*n + (i+1)*n + j;
          v[3] = k*n*n + (i+1)*n + j+1;

          v[4] = (k+1)*n*n + i*n + j;
          v[5] = (k+1)*n*n + i*n + j+1;
          v[6] = (k+1)*n*n + (i+1)*n + j;
          v[7] = (k+1)*n*n + (i+1)*n + j+1;

          factory.insertElement(GeometryType(GeometryType::cube,dim), v);
        }
      }
    }
  }
  else if (dim == 2) {
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
  }
};

template <class GridType>
void createGrid(GridType &grid, int n)
{
  typedef Dune::GridFactory<GridType> GridFactory;

  GridFactory factory(&grid);

  insertVertices<GridType>(n, factory);
  insertElements<GridType>(n, factory);

  factory.createGrid();
}

template <class GridType, int commCodim>
void testCommunication(GridType &grid)
{
  std::cout << "Testing communication for codim " << commCodim << " entities\n";

  typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
      LayoutWrapper<commCodim>::template Layout> LeafMapperType;

  typedef LevelMultipleCodimMultipleGeomTypeMapper<GridType,
      LayoutWrapper<commCodim>::template Layout> LevelMapperType;

  LevelMapperType levelMapper(grid, 0);
  LeafMapperType leafMapper(grid);

  std::cout << "Level index set has " << levelMapper.size() << " codim " << commCodim << " entities\n";
  std::cout << "Leaf index set has " << leafMapper.size() << " codim " << commCodim << " entities\n";

  assert(levelMapper.size() == leafMapper.size());

  const int dim = GridType::dimension;

  const typename GridType::LeafGridView &leafGridView = grid.leafView();

  // create the user data arrays
  typedef std::vector<Dune::FieldVector<double, 1> > UserDataType;
  UserDataType userDataSend(leafGridView.size(commCodim), 0.0);
  UserDataType userDataReceive(leafGridView.size(commCodim), 0.0);
  UserDataType entityIndex(leafGridView.size(commCodim), -1e10);
  UserDataType partitionType(leafGridView.size(commCodim), -1e10);

  // write the partition type of each entity into the corrosponding
  // result array
  typename GridType::LeafGridView::template Codim<commCodim>::Iterator
  it = leafGridView.template begin<commCodim>();
  const typename GridType::LeafGridView::template Codim<commCodim>::Iterator
  &endIt = leafGridView.template end<commCodim>();
  for (; it != endIt; ++it) {
    entityIndex[leafMapper.map(*it)] = leafMapper.map(*it);
    partitionType[leafMapper.map(*it)] = it->partitionType();
  };

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
  // Write results to disk
  //////////////////////////////////////////////////////
  std::cout << "Writing data to disk\n";
  Dune::VTKWriter<typename GridType::LeafGridView> writer(leafGridView);
  if (commCodim == 0) {
    writer.addCellData(userDataSend, "Send");
    writer.addCellData(userDataReceive, "Receive");
    writer.addCellData(partitionType, "Partition Type");
    writer.addCellData(entityIndex, "Entity Index");
  }
  else if (commCodim == dim) {
    writer.addVertexData(userDataSend, "Send");
    writer.addVertexData(userDataReceive, "Receive");
    writer.addVertexData(partitionType, "Partition Type");
    writer.addVertexData(entityIndex, "Entity Index");
  }

  char fileName[1024];
  sprintf(fileName, "test-parallel-ug-dim=%d-commCodim=%d", dim, commCodim);
  writer.write(fileName,
               Dune::VTKOptions::ascii);
  std::cout << "Done writing data to disk\n";
};

template <int dim>
void testParallelUG()
{
  std::cout << "Testing parallel UGGrid for " << dim << "D\n";

  ////////////////////////////////////////////////////////////
  // Create uniform grid on rank 0
  ////////////////////////////////////////////////////////////
  int n = 11;
  typedef UGGrid<dim> GridType;
  GridType grid;
  createGrid(grid, n);


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


  std::cout << "Process " << grid.comm().rank() + 1
            << " has " << grid.size(0)
            << " elements and " << grid.size(dim)
            << " nodes.\n";

  // make sure each process has some elements
  assert(grid.size(0) > 0);
  // make sure each process has some nodes
  assert(grid.size(dim) > 0);

  //////////////////////////////////////////////////////
  // Test level and leaf mappers/gridViews
  //////////////////////////////////////////////////////
  typename GridType::LevelGridView level0GridView = grid.levelView(0);
  typename GridType::LeafGridView leafGridView = grid.leafView();

  std::cout << "LevelGridView for level 0 has " << level0GridView.size(0)
            << " elements "  << level0GridView.size(dim - 1)
            << " edges and " << level0GridView.size(dim)
            << " nodes.\n";
  std::cout << "LeafGridView has " << leafGridView.size(0)
            << " elements, " << leafGridView.size(dim - 1)
            << " edges and " << leafGridView.size(dim)
            << " nodes.\n";

  // iterate over all elements and element intersections
  checkConsistency(level0GridView);
  checkConsistency(leafGridView);

  //////////////////////////////////////////////////////
  // Test element and node communication
  //////////////////////////////////////////////////////
  testCommunication<GridType, 0>(grid);
  testCommunication<GridType, dim>(grid);

  // //////////////////////////////////////////////////
  //   Refine globally and test again
  // //////////////////////////////////////////////////

  grid.globalRefine(1);

  for (int i=0; i<=grid.maxLevel(); i++)
    checkConsistency(grid.levelView(i));

  checkConsistency(grid.leafView());

  testCommunication<GridType, 0>(grid);
  testCommunication<GridType, dim>(grid);

};

int main (int argc , char **argv) try
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper &mpiHelper =
    Dune::MPIHelper::instance(argc, argv);

  std::cout << "This is process "
            << mpiHelper.rank() + 1
            << " of "
            << mpiHelper.size()
            << ", PID "
            << getpid()
            << " .\n";

  // test 2D grid
  testParallelUG<2>();

  // test 3D grid
  testParallelUG<3>();

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
