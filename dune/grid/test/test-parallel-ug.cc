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
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/geometry/referenceelements.hh>

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
void checkIntersections(const GridView &gv)
{
  typedef typename GridView::template Codim<0>::Entity Element;
  typedef typename GridView::template Codim<0>::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  Iterator it = gv.template begin<0>();
  const Iterator &endIt = gv.template end<0>();
  // check the intersections
  for (; it != endIt; ++it) {
    if (it->partitionType() == Dune::GhostEntity)
      continue;

    IntersectionIterator isIt           = gv.ibegin(*it);
    const IntersectionIterator &isEndIt = gv.iend(*it);
    int n = 0;
    for (; isIt != isEndIt; ++isIt) {
      isIt->boundary();
      isIt->inside();
      if (isIt->neighbor()) {
        isIt->outside();
      }

      ++ n;
    }

    if (n != it->template count<1>())
    {

      DUNE_THROW(Dune::InvalidStateException,
                 "Number of faces for non-ghost cell incorrect. Is "
                 << n << " but should be " << it->template count<1>());
    }
  }
}

// some consistency tests for the mappers
template <int codim, class GridView>
void checkMappers(const GridView &gridView)
{
  typename GridView::template Codim<codim>::Iterator
  it = gridView.template begin<codim>();
  const typename GridView::template Codim<codim>::Iterator
  &endIt = gridView.template end<codim>();

  // make sure the number of entities is the same for
  // gridView.size() and iterating through the grid
  int numEntities = 0;
  for (; it != endIt; ++it)
    ++ numEntities;
  if (numEntities != gridView.size(codim)) {
    DUNE_THROW(InvalidStateException,
               "Number of codim " << codim
                                  << " entities is inconsistent (iterator: " << numEntities
                                  << " grid view: " << gridView.size(codim) << ")");
  };

  typedef Dune::MultipleCodimMultipleGeomTypeMapper
  <GridView, LayoutWrapper<codim>::template Layout> MapperType;
  MapperType mapper(gridView);

  // make sure no entity has two indices
  std::vector<int> indices(numEntities, -100);
  it = gridView.template begin<codim>();
  for (; it != endIt; ++it) {
    int i = mapper.map(*it);
    if (i < 0 || i >= numEntities) {
      DUNE_THROW(InvalidStateException,
                 "Mapper returns an invalid index for codim " << codim
                                                              << " entities: " << i
                                                              << " should be in range [0,"
                                                              << numEntities -1 <<  "]");
    };
    if (indices[i] >= 0) {
      DUNE_THROW(InvalidStateException,
                 "Mapper returns index " << i << " twice for codim " << codim
                                         << " entities.");
    }
    indices[i] = i;
  };
};

// specializations for non-implementet cases
template <int dim, int codim, class GridView>
struct checkMappersWrapper
{
  static void check(const GridView &gv)
  { return checkMappers<codim>(gv); }
};

// codim 3 entities for 2d grids
template <class GridView>
struct checkMappersWrapper<2, 3, GridView>
{ static void check(const GridView &gv) { } };

// codim 1 entities for 3d grids
template <class GridView>
struct checkMappersWrapper<3, 1, GridView>
{ static void check(const GridView &gv) { } };

// codim 2 entities for 3d grids
template <class GridView>
struct checkMappersWrapper<3, 2, GridView>
{ static void check(const GridView &gv) { } };

// codim 1 entities for 2d grids
template <class GridView>
struct checkMappersWrapper<2, 1, GridView>
{ static void check(const GridView &gv) { } };


template <class GridView, int commCodim>
void testCommunication(const GridView &gridView, bool isLeaf, bool printVTK=false)
{
  std::cout << "Testing communication for codim " << commCodim << " entities\n";

  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, LayoutWrapper<commCodim>::template Layout>
  MapperType;
  MapperType mapper(gridView);

  std::cout << "Index set has " << mapper.size() << " codim " << commCodim << " entities\n";

  const int dim = GridView::dimension;

  // create the user data arrays
  typedef std::vector<Dune::FieldVector<double, 1> > UserDataType;
  UserDataType userDataSend(gridView.size(commCodim), 0.0);
  UserDataType userDataReceive(gridView.size(commCodim), 0.0);
  UserDataType entityIndex(gridView.size(commCodim), -1e10);
  UserDataType partitionType(gridView.size(commCodim), -1e10);

  // write the partition type of each entity into the corrosponding
  // result array
  typename GridView::template Codim<commCodim>::Iterator
  it = gridView.template begin<commCodim>();
  const typename GridView::template Codim<commCodim>::Iterator
  &endIt = gridView.template end<commCodim>();
  for (; it != endIt; ++it) {
    entityIndex[mapper.map(*it)]   = mapper.map(*it);
    partitionType[mapper.map(*it)] = it->partitionType();
  };

  // initialize data handle (marks the nodes where some data was
  // send or received)
  typedef DataExchange<MapperType, commCodim> MyDataHandle;
  MyDataHandle datahandle(mapper,
                          userDataSend,
                          userDataReceive);

  // communicate the entities at the interior border to all other
  // processes
  gridView.communicate(datahandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

  //////////////////////////////////////////////////////
  // Write results to disk
  //////////////////////////////////////////////////////
  if (printVTK)
  {
    std::cout << "Writing data to disk\n";
    Dune::VTKWriter<GridView> writer(gridView);
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
    if (isLeaf)
      strcat(fileName, "-leaf");
    else
    {
      char levelName[10];
      sprintf(levelName, "-level=%d", (gridView.template begin<commCodim>())->level());
      strcat(fileName, levelName);
    }
    writer.write(fileName, Dune::VTK::ascii);
    std::cout << "Done writing data to disk\n";
  }
}

//! edge communication
template <int dim>
class EdgeCommunication
{
public:
  template <class GridView>
  static void test(const GridView &gridView)
  {
    const int dim = GridView::dimension;
    const int commCodim = dim - 1;

    std::cout << "Testing communication for codim " << commCodim << " entities\n";

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, LayoutWrapper<commCodim>::template Layout>
    MapperType;
    MapperType mapper(gridView);

    std::cout << gridView.comm().rank()+1 << ": Index set has " << mapper.size() << " codim " << commCodim << " entities\n";

    // create the user data arrays
    typedef std::vector<Dune::FieldVector<double, 1> > UserDataType;
    UserDataType userDataSend(gridView.size(commCodim), 0.0);
    UserDataType userDataReceive(gridView.size(commCodim), 0.0);
    UserDataType entityIndex(gridView.size(commCodim), -1e10);
    UserDataType partitionType(gridView.size(commCodim), -1e10);

    // write the partition type of each entity into the corresponding
    // result array
    typename GridView::template Codim<0>::Iterator
    it = gridView.template begin<0>();
    const typename GridView::template Codim<0>::Iterator
    &endIt = gridView.template end<0>();
    for (; it != endIt; ++it) {
      int numberOfEdges = it->template count<commCodim>();
      for (int k = 0; k < numberOfEdges; k++)
      {
        typedef typename GridView::template Codim<0>::Entity Element;
        typedef typename Element::template Codim<commCodim>::EntityPointer EdgePointer;
        const EdgePointer edgePointer(it->template subEntity<commCodim>(k));
        entityIndex[mapper.map(*edgePointer)]   = mapper.map(*edgePointer);
        partitionType[mapper.map(*edgePointer)] = edgePointer->partitionType();

        if (edgePointer->partitionType() == Dune::BorderEntity)
        {
          typedef typename GridView::template Codim<0>::Entity Element;
          const typename Element::Geometry& geometry = it->geometry();
          Dune::GeometryType gt = geometry.type();

          const typename Dune::GenericReferenceElementContainer<double, dim>::value_type&
          referenceElement = Dune::GenericReferenceElements<double, dim>::general(gt);
          const Dune::FieldVector<double, dim>&
          edgeGlobal = geometry.global(referenceElement.position(k, dim-1));
          std::cout << gridView.comm().rank()+1 << ": border edge "
                    << mapper.map(*edgePointer) << " (" << edgeGlobal
                    << ")" << std::endl;
        }
      }
    }

    // initialize data handle (marks the nodes where some data was
    // send or received)
    typedef DataExchange<MapperType, commCodim> MyDataHandle;
    MyDataHandle datahandle(mapper,
                            userDataSend,
                            userDataReceive);

    // communicate the entities at the interior border to all other
    // processes
    gridView.communicate(datahandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
  }
};


template <int dim>
void testParallelUG()
{
  std::cout << "Testing parallel UGGrid for " << dim << "D\n";

  ////////////////////////////////////////////////////////////
  // Create uniform grid on rank 0
  ////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  StructuredGridFactory<GridType> structuredGridFactory;

  Dune::FieldVector<double,dim> lowerLeft(0);
  Dune::FieldVector<double,dim> upperRight(1);
  Dune::array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 4);
  shared_ptr<GridType> grid = structuredGridFactory.createCubeGrid(lowerLeft, upperRight, elements);

  //////////////////////////////////////////////////////
  // Distribute the grid
  //////////////////////////////////////////////////////
  /*    grid.loadBalance(0,   // strategy
                       0,   // minlevel
                       2,   // depth
                       32,   // maxlevel
                       1);   // minelement
   */
  grid->loadBalance();


  std::cout << "Process " << grid->comm().rank() + 1
            << " has " << grid->size(0)
            << " elements and " << grid->size(dim)
            << " nodes.\n";

  // make sure each process has some elements
  assert(grid->size(0) > 0);
  // make sure each process has some nodes
  assert(grid->size(dim) > 0);

  //////////////////////////////////////////////////////
  // Test level and leaf mappers/gridViews
  //////////////////////////////////////////////////////

  typedef typename GridType::LevelGridView LevelGV;
  typedef typename GridType::LeafGridView LeafGV;
  const LeafGV&  leafGridView = grid->leafView();
  const LevelGV& level0GridView = grid->levelView(0);

  std::cout << "LevelGridView for level 0 has " << level0GridView.size(0)
            << " elements "  << level0GridView.size(dim - 1)
            << " edges and " << level0GridView.size(dim)
            << " nodes.\n";
  std::cout << "LeafGridView has " << leafGridView.size(0)
            << " elements, " << leafGridView.size(dim - 1)
            << " edges and " << leafGridView.size(dim)
            << " nodes.\n";

  // some consistency checks for the mappers and the intersections
  checkIntersections(level0GridView);
  checkMappersWrapper<dim, 0, LevelGV>::check(level0GridView);
  checkMappersWrapper<dim, 1, LevelGV>::check(level0GridView);
  checkMappersWrapper<dim, 2, LevelGV>::check(level0GridView);
  checkMappersWrapper<dim, 3, LevelGV>::check(level0GridView);

  checkIntersections(leafGridView);
  checkMappersWrapper<dim, 0, LeafGV>::check(leafGridView);
  checkMappersWrapper<dim, 1, LeafGV>::check(leafGridView);
  checkMappersWrapper<dim, 2, LeafGV>::check(leafGridView);
  checkMappersWrapper<dim, 3, LeafGV>::check(leafGridView);

  // Test element and node communication on level view
  testCommunication<typename GridType::LevelGridView, 0>(level0GridView, false);
  testCommunication<typename GridType::LevelGridView, dim>(level0GridView, false);
  EdgeCommunication<dim>::template test<typename GridType::LevelGridView>(level0GridView);

  // Test element and node communication on leaf view
  testCommunication<typename GridType::LeafGridView, 0>(leafGridView, true);
  testCommunication<typename GridType::LeafGridView, dim>(leafGridView, true);
  EdgeCommunication<dim>::template test<typename GridType::LeafGridView>(leafGridView);

  ////////////////////////////////////////////////////
  //  Refine globally and test again
  ////////////////////////////////////////////////////

  std::cout << "Global refinement\n";
  grid->globalRefine(1);

  for (int i=0; i<=grid->maxLevel(); i++) {
    checkIntersections(grid->levelView(i));
    checkMappersWrapper<dim, 0, LevelGV>::check(grid->levelView(i));
    checkMappersWrapper<dim, 1, LevelGV>::check(grid->levelView(i));
    checkMappersWrapper<dim, 2, LevelGV>::check(grid->levelView(i));
    checkMappersWrapper<dim, 3, LevelGV>::check(grid->levelView(i));
  }

  checkIntersections(grid->leafView());
  checkMappersWrapper<dim, 0, LeafGV>::check(grid->leafView());
  checkMappersWrapper<dim, 1, LeafGV>::check(grid->leafView());
  checkMappersWrapper<dim, 2, LeafGV>::check(grid->leafView());
  checkMappersWrapper<dim, 3, LeafGV>::check(grid->leafView());

  for (int i=0; i<=grid->maxLevel(); i++)
  {
    testCommunication<typename GridType::LevelGridView, 0>(grid->levelView(i), false);
    testCommunication<typename GridType::LevelGridView, dim>(grid->levelView(i), false);
    EdgeCommunication<dim>::template test<typename GridType::LevelGridView>(grid->levelView(i));
  }
  testCommunication<typename GridType::LeafGridView, 0>(grid->leafView(), true);
  testCommunication<typename GridType::LeafGridView, dim>(grid->leafView(), true);
  EdgeCommunication<dim>::template test<typename GridType::LeafGridView>(grid->leafView());

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
