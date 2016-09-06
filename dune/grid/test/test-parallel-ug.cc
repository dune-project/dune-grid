// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Test parallel interface if a parallel UG is used

#include <config.h>

#include <unistd.h>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

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
  bool fixedSize (int dim, int codim) const
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
    DataType x = mapper_.index(e);

    userDataSend_[mapper_.index(e)][0] = x;
    std::cout << "Process "
              << Dune::MPIHelper::getCollectiveCommunication().rank()+1
              << " sends for entity "
              << mapper_.index(e)
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

    userDataReceive_[mapper_.index(e)][0] = x;
    std::cout << "Process "
              << Dune::MPIHelper::getCollectiveCommunication().rank()+1
              << " received for entity "
              << mapper_.index(e)
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
  for (const auto& element : elements(gv)) {
    if (element.partitionType() == Dune::GhostEntity)
      continue;

    unsigned int n = 0;
    for (const auto& intersection : intersections(gv, element)) {
      intersection.boundary();
      intersection.inside();
      if (intersection.neighbor()) {
        intersection.outside();
      }

      ++ n;
    }

    if (n != element.subEntities(1))
    {

      DUNE_THROW(Dune::InvalidStateException,
                 "Number of faces for non-ghost cell incorrect. Is "
                 << n << " but should be " << element.subEntities(1));
    }
  }
}

// some consistency tests for the mappers
template <int codim, class GridView>
void checkMappers(const GridView &gridView)
{
  // make sure the number of entities is the same for
  // gridView.size() and iterating through the grid
  int numEntities = 0;
  for (DUNE_UNUSED const auto& entity : entities(gridView, Dune::Codim<codim>()))
    ++ numEntities;
  if (numEntities != gridView.size(codim)) {
    DUNE_THROW(InvalidStateException,
               gridView.comm().rank() + 1
               << ": Number of codim " << codim
               << " entities is inconsistent (iterator: " << numEntities
               << " grid view: " << gridView.size(codim) << ")");
  }

  typedef Dune::MultipleCodimMultipleGeomTypeMapper
  <GridView, LayoutWrapper<codim>::template Layout> MapperType;
  MapperType mapper(gridView);

  // make sure no entity has two indices
  std::vector<int> indices(numEntities, -100);
  for (const auto& entity : entities(gridView, Dune::Codim<codim>())) {
    int i = mapper.index(entity);
    if (i < 0 || i >= numEntities) {
      DUNE_THROW(InvalidStateException,
                 gridView.comm().rank() + 1
                 << ": Mapper returns an invalid index for codim " << codim
                 << " entities: " << i
                 << " should be in range [0,"
                 << numEntities -1 <<  "]");
    }
    if (indices[i] >= 0) {
      DUNE_THROW(InvalidStateException,
                 gridView.comm().rank() + 1
                 << ": Mapper returns index " << i << " twice for codim " << codim
                 << " entities.");
    }
    indices[i] = i;
  }
}

// specializations for non-implemented cases
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
  std::cout << gridView.comm().rank() + 1
            << ": Testing communication for codim " << commCodim << " entities\n";

  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, LayoutWrapper<commCodim>::template Layout>
  MapperType;
  MapperType mapper(gridView);

  std::cout << gridView.comm().rank() + 1
            << ": Index set has " << mapper.size() << " codim " << commCodim << " entities\n";

  const int dim = GridView::dimension;

  // create the user data arrays
  typedef std::vector<Dune::FieldVector<double, 1> > UserDataType;
  UserDataType userDataSend(gridView.size(commCodim), 0.0);
  UserDataType userDataReceive(gridView.size(commCodim), 0.0);
  UserDataType entityIndex(gridView.size(commCodim), -1e10);
  UserDataType partitionType(gridView.size(commCodim), -1e10);

  // write the partition type of each entity into the corrosponding
  // result array
  for (const auto& entity : entities(gridView, Dune::Codim<commCodim>())) {
    entityIndex[mapper.index(entity)]   = mapper.index(entity);
    partitionType[mapper.index(entity)] = entity.partitionType();
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

//! edge and face communication
template <class GridView, int commCodim>
class EdgeAndFaceCommunication
{
public:
  static void test(const GridView &gridView)
  {
    const int dim = GridView::dimension;

    std::cout << gridView.comm().rank() + 1
              << ": Testing communication for codim " << commCodim << " entities\n";

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
    for (const auto& element : elements(gridView)) {
      int numberOfSubEntities = element.subEntities(commCodim);
      for (int k = 0; k < numberOfSubEntities; k++)
      {
        typedef typename GridView::template Codim<0>::Entity Element;
        typedef typename Element::template Codim<commCodim>::Entity Entity;
        const Entity entity(element.template subEntity<commCodim>(k));
        entityIndex[mapper.index(entity)]   = mapper.index(entity);
        partitionType[mapper.index(entity)] = entity.partitionType();

        if (entity.partitionType() == Dune::BorderEntity)
        {
          const auto geometry = element.geometry();
          const auto gt = geometry.type();

          const auto& referenceElement = Dune::ReferenceElements<double, dim>::general(gt);
          const auto entityGlobal = geometry.global(referenceElement.position(k, commCodim));
          std::cout << gridView.comm().rank()+1 << ": border codim "
                    << commCodim << " entity "
                    << mapper.index(entity) << " (" << entityGlobal
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

    std::cout << gridView.comm().rank() + 1
              << ": Finished testing communication for codim " << commCodim << " entities\n";
  }
};

class LoadBalance
{
  template<class Grid, class Vector, int commCodim>
  class LBDataHandle
    : public Dune::CommDataHandleIF<LBDataHandle<Grid, Vector, commCodim>,
          typename Vector::value_type>
  {
  public:
    typedef typename Vector::value_type DataType;
    typedef Dune::CommDataHandleIF<LBDataHandle<Grid, Vector, commCodim>,
        DataType> ParentType;

    bool contains (int dim, int codim) const
    { return (codim == commCodim); }

    bool fixedSize (int dim, int codim) const
    { return true; }

    template<class Entity>
    size_t size (Entity& entity) const
    {
      return 1;
    }

    template<class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& entity) const
    {
      int index = grid_.leafGridView().indexSet().index(entity);
      buff.write(dataVector_[index]);
    }

    template<class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& entity, size_t n)
    {
      if (dataVector_.size() != grid_.leafGridView().size(commCodim))
        dataVector_.resize(grid_.leafGridView().size(commCodim));

      int index = grid_.leafGridView().indexSet().index(entity);
      buff.read(dataVector_[index]);
    }

    LBDataHandle (Grid& grid, Vector& dataVector)
      : grid_(grid), dataVector_(dataVector)
    {}

  private:
    Grid& grid_;
    Vector& dataVector_;
  };

public:
  template <class Grid>
  static void test(Grid& grid)
  {
    const int dim = Grid::dimension;
    const int commCodim = dim;
    typedef typename Grid::ctype ctype;

    // define the vector containing the data to be balanced
    typedef Dune::FieldVector<ctype, dim> Position;
    std::vector<Position> dataVector(grid.leafGridView().size(commCodim));

    // fill the data vector
    const auto& gv = grid.leafGridView();
    for (const auto& entity : entities(gv, Dune::Codim<commCodim>(),
                                       Dune::Partitions::interiorBorder)) {
      int index = gv.indexSet().index(entity);

      // assign the position of the entity to the entry in the vector
      dataVector[index] = entity.geometry().center();
    }

    // balance the grid and the data
    LBDataHandle<Grid, std::vector<Position>, commCodim> dataHandle(grid, dataVector);
    grid.loadBalance(dataHandle);

    // check for correctness
    for (const auto& entity : entities(gv, Dune::Codim<commCodim>(),
                                       Dune::Partitions::interiorBorder)) {
      int index = gv.indexSet().index(entity);

      const auto position = entity.geometry().center();

      // compare the position with the balanced data
      for (int k = 0; k < dim; k++)
      {
        if (Dune::FloatCmp::ne(dataVector[index][k], position[k]))
        {
          DUNE_THROW(Dune::ParallelError,
                     gv.comm().rank() << ": position " << position
                                      << " does not coincide with communicated data "
                                      << dataVector[index]);
        }
      }
    }

    std::cout << gv.comm().rank()
              << ": load balancing with data was successful." << std::endl;
  }
};

template <int dim>
void testParallelUG(bool simplexGrid, bool localRefinement, int refinementDim, bool refineUpperPart)
{
  std::cout << "Testing parallel UGGrid for " << dim << "D\n";

  ////////////////////////////////////////////////////////////
  // Create uniform grid on rank 0
  ////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  StructuredGridFactory<GridType> structuredGridFactory;

  Dune::FieldVector<double,dim> lowerLeft(0);
  Dune::FieldVector<double,dim> upperRight(1);
  std::array<unsigned int, dim> numElements;
  std::fill(numElements.begin(), numElements.end(), 4);
  std::shared_ptr<GridType> grid;
  if (simplexGrid)
    grid = structuredGridFactory.createSimplexGrid(lowerLeft, upperRight, numElements);
  else
    grid = structuredGridFactory.createCubeGrid(lowerLeft, upperRight, numElements);

  //////////////////////////////////////////////////////
  // Distribute the grid
  //////////////////////////////////////////////////////
  LoadBalance::test(*grid);

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
  const LeafGV&  leafGridView = grid->leafGridView();
  const LevelGV& level0GridView = grid->levelGridView(0);

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
  testCommunication<LevelGV, 0>(level0GridView, false);
  testCommunication<LevelGV, dim>(level0GridView, false);
  EdgeAndFaceCommunication<LevelGV, dim-1>::test(level0GridView);
  if (dim == 3)
    EdgeAndFaceCommunication<LevelGV, 1>::test(level0GridView);

  // Test element and node communication on leaf view
  testCommunication<LeafGV, 0>(leafGridView, true);
  testCommunication<LeafGV, dim>(leafGridView, true);
  EdgeAndFaceCommunication<LeafGV, dim-1>::test(leafGridView);
  if (dim == 3)
    EdgeAndFaceCommunication<LeafGV, 1>::test(leafGridView);

  ////////////////////////////////////////////////////
  //  Refine globally and test again
  ////////////////////////////////////////////////////

  if (!localRefinement)
  {
    std::cout << "Global refinement\n";
    grid->globalRefine(1);
  }
  else
  {
    std::cout << "Local refinement\n";

    // mark all elements with x-coordinate < 0.5 for refinement
    for (const auto& element : elements(grid->leafGridView())) {
      int nRefine = 1;
      bool refine;
      if (refineUpperPart) {
        refine = (element.geometry().center()[refinementDim] > 0.5);
      }
      else {
        refine = (element.geometry().center()[refinementDim] < 0.5);
      }
      if (refine)
        grid->mark(nRefine, element);
    }

    // adapt the grid
    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();

    // write the adapted grid to VTK
    Dune::VTKWriter<LeafGV> writer(grid->leafGridView());
    char fileName[1024];
    sprintf(fileName, "adapted-grid-dim=%d", dim);
    writer.write(fileName, Dune::VTK::ascii);
  }

  for (int i=0; i<=grid->maxLevel(); i++) {
    checkIntersections(grid->levelGridView(i));
    checkMappersWrapper<dim, 0, LevelGV>::check(grid->levelGridView(i));
    checkMappersWrapper<dim, 1, LevelGV>::check(grid->levelGridView(i));
    checkMappersWrapper<dim, 2, LevelGV>::check(grid->levelGridView(i));
    checkMappersWrapper<dim, 3, LevelGV>::check(grid->levelGridView(i));
  }

  checkIntersections(grid->leafGridView());
  checkMappersWrapper<dim, 0, LeafGV>::check(grid->leafGridView());
  checkMappersWrapper<dim, 1, LeafGV>::check(grid->leafGridView());
  checkMappersWrapper<dim, 2, LeafGV>::check(grid->leafGridView());
  checkMappersWrapper<dim, 3, LeafGV>::check(grid->leafGridView());

  for (int i=0; i<=grid->maxLevel(); i++)
  {
    testCommunication<LevelGV, 0>(grid->levelGridView(i), false);
    testCommunication<LevelGV, dim>(grid->levelGridView(i), false);
    EdgeAndFaceCommunication<LevelGV, dim-1>::test(grid->levelGridView(i));
    if (dim == 3)
      EdgeAndFaceCommunication<LevelGV, 1>::test(grid->levelGridView(i));
  }
  testCommunication<LeafGV, 0>(grid->leafGridView(), true);
  testCommunication<LeafGV, dim>(grid->leafGridView(), true);
  EdgeAndFaceCommunication<LeafGV, dim-1>::test(grid->leafGridView());
  if (dim == 3)
    EdgeAndFaceCommunication<LeafGV, 1>::test(grid->leafGridView());

}

int main (int argc , char **argv) try
{
  // initialize MPI, finalize is done automatically on exit
  auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

  std::cout << "This is process "
            << mpiHelper.rank() + 1
            << " of "
            << mpiHelper.size()
            << ", PID "
            << getpid()
            << " .\n";

  /*
   * Test 2D and 3D grids,
   * on structured cube and simplex grids,
   * with global and local refinement.
   */
  for (const bool simplexGrid : {false, true}) {
    for (const bool localRefinement : {false, true}) {
      for (const bool refineUpperPart : {false, true}) {
        for (const int refinementDim : {0,1})
          testParallelUG<2>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
        for (const int refinementDim : {0,1,2})
          testParallelUG<3>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
      }
    }
  }

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
