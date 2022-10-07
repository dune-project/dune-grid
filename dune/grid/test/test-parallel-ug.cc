// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Test parallel interface if a parallel UG is used

#include <config.h>

#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <bitset>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/stdstreams.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;

// A DataHandle class to exchange entries of a vector.
template<class MapperT>
class DataExchange
  : public Dune::CommDataHandleIF<DataExchange<MapperT>,
        double>
{
public:
  //! export type of data for message buffer
  typedef double DataType;
  typedef std::vector<Dune::FieldVector<double,1> > UserDataType;

  //! constructor
  DataExchange(const MapperT &mapper,
               std::bitset<MapperT::GridView::dimension+1> communicationCodims,
               std::vector<std::size_t>& gatherCounter,
               std::vector<std::size_t>& scatterCounter,
               UserDataType &userDataSend,
               UserDataType &userDataReceive)
    : mapper_(mapper),
      communicationCodims_(communicationCodims),
      gatherCounter_(gatherCounter),
      scatterCounter_(scatterCounter),
      userDataSend_(userDataSend),
      userDataReceive_(userDataReceive)
  {}


  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return communicationCodims_[codim];
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
    gatherCounter_[mapper_.index(e)]++;

    userDataSend_[mapper_.index(e)][0] = x;
    dverb << "Process "
              << Dune::MPIHelper::getCommunication().rank()+1
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
    scatterCounter_[mapper_.index(e)]++;
    dverb << "Process "
              << Dune::MPIHelper::getCommunication().rank()+1
              << " received for entity "
              << mapper_.index(e)
              << ": "
              << std::setprecision(20)
              << x << "\n";
  }

private:
  const MapperT &mapper_;
  const std::bitset<MapperT::GridView::dimension+1> communicationCodims_;
  std::vector<std::size_t> &gatherCounter_;
  std::vector<std::size_t> &scatterCounter_;
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
  for ([[maybe_unused]] const auto& entity : entities(gridView, Dune::Codim<codim>()))
    ++ numEntities;
  if (numEntities != gridView.size(codim)) {
    DUNE_THROW(InvalidStateException,
               gridView.comm().rank() + 1
               << ": Number of codim " << codim
               << " entities is inconsistent (iterator: " << numEntities
               << " grid view: " << gridView.size(codim) << ")");
  }

  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> MapperType;
  MapperType mapper(gridView, mcmgLayout(Codim<codim>{}));

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


template <class GridView>
void testCommunication(const GridView &gridView,
                       std::bitset<GridView::dimension+1> communicationCodims,
                       InterfaceType communicationInterface,
                       CommunicationDirection communicationDirection,
                       const std::set<PartitionType>& sendingPartitions,
                       const std::set<PartitionType>& receivingPartitions)
{
  const int dim = GridView::dimension;

  dverb << gridView.comm().rank() + 1
            << ": Testing communication for the following codimensions " << communicationCodims << std::endl;

  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> MapperType;
  auto mcmgLayout = [&communicationCodims](GeometryType gt, int dimgrid)
    {
      return communicationCodims[dim - gt.dim()];
    };
  MapperType mapper(gridView, mcmgLayout);

  // create the user data arrays
  using UserDataType = std::vector<Dune::FieldVector<double, 1> >;

  UserDataType userDataSend(mapper.size(), 0.0);
  UserDataType userDataReceive(mapper.size(), 0.0);
  // For each entity, count how many times 'gather' and 'scatter' have been called
  std::vector<std::size_t> gatherCounter(mapper.size(), 0);
  std::vector<std::size_t> scatterCounter(mapper.size(), 0);

  // write the partition type of each entity into the corresponding
  // result array
  for (const auto& element : elements(gridView))
  {
    Hybrid::forEach(std::make_index_sequence< dim+1 >{},
      [&](auto codim){
      if (communicationCodims[codim])
      {
        auto numberOfSubEntities = element.subEntities(codim);
        for (std::size_t k = 0; k < numberOfSubEntities; k++)
        {
          const auto entity(element.template subEntity<codim>(k));

          if (entity.partitionType() == Dune::BorderEntity)
          {
            const auto geometry = element.geometry();

            auto referenceElement = Dune::referenceElement<double, dim>(element.type());
            const auto entityGlobal = geometry.global(referenceElement.position(k, codim));
            dverb << gridView.comm().rank()+1 << ": border codim "
                  << codim << " entity "
                  << mapper.index(entity) << " (" << entityGlobal
                  << ")" << std::endl;
          }
        }
      }
    });
  }

  // initialize data handle (marks the nodes where some data was
  // send or received)
  DataExchange<MapperType> datahandle(mapper,
                                      communicationCodims,
                                      gatherCounter,
                                      scatterCounter,
                                      userDataSend,
                                      userDataReceive);

  // communicate the entities at the interior border to all other
  // processes
  gridView.communicate(datahandle, communicationInterface, communicationDirection);

  // write the partition type of each entity into the corresponding
  // result array
  for (const auto& element : elements(gridView))
  {
    Hybrid::forEach(std::make_index_sequence< dim+1 >{},
      [&](auto codim)
    {
      // TODO: Also check whether no communication happens when no communication
      // is requested.
      if (communicationCodims[codim])
      {
        auto numberOfSubEntities = element.subEntities(codim);
        for (std::size_t k = 0; k < numberOfSubEntities; k++)
        {
          const auto entity(element.template subEntity<codim>(k));

          auto partitionTypes = entity.impl().partitionTypes();

          // Check whether 'gather' has been called the appropriate number of times
          std::size_t expectedNumberOfGatherCalls = 0;

          for (const auto& pType : partitionTypes)
            if (pType.first != gridView.comm().rank()
                && (sendingPartitions.find(entity.partitionType())!=sendingPartitions.end())
                && (receivingPartitions.find(pType.second))!=receivingPartitions.end())
              expectedNumberOfGatherCalls++;

          if (gatherCounter[mapper.index(entity)] != expectedNumberOfGatherCalls)
          {
            std::cerr << gridView.comm().rank() << ": UGGrid did not call 'gather' "
                      << expectedNumberOfGatherCalls << " times, but "
                      << gatherCounter[mapper.index(entity)] << " times on an entity!" << std::endl;
            std::cerr << gridView.comm().rank() << ": Problematic entity: codim = " << codim
                      << ",  partitionType = " << entity.partitionType()
                      << ",  center = " << entity.geometry().center()
                      << std::endl;
            std::abort();
          }

          // Check whether 'scatter' has been called the appropriate number of times
          std::size_t expectedNumberOfScatterCalls = 0;

          for (const auto& pType : partitionTypes)
            if (pType.first != gridView.comm().rank()
                && (receivingPartitions.find(entity.partitionType()) != receivingPartitions.end())
                && (sendingPartitions.find(pType.second) != sendingPartitions.end()))
              expectedNumberOfScatterCalls++;

          if (scatterCounter[mapper.index(entity)] != expectedNumberOfScatterCalls)
          {
            std::cerr << gridView.comm().rank() << ": UGGrid did not call 'scatter' "
                      << expectedNumberOfScatterCalls << " times, but "
                      << scatterCounter[mapper.index(entity)] << " times on an entity!" << std::endl;
            std::cerr << gridView.comm().rank() << ": Problematic entity: codim = " << codim
                      << ",  partitionType = " << entity.partitionType()
                      << ",  center = " << entity.geometry().center()
                      << std::endl;
            std::abort();
          }
        }
      }
    });
  }
}

template<typename Grid>
class LoadBalance
{
  const static int dimension = Grid::dimension;
  using ctype = typename Grid::ctype;
  using Position = Dune::FieldVector<ctype, dimension>;
  using GridView = typename Grid::LeafGridView;
  using IdSet = typename Grid::LocalIdSet;
  using Data = std::map<typename IdSet::IdType, Position>;

  using Codims = std::bitset<dimension+1>;

  class LBDataHandle
    : public Dune::CommDataHandleIF<LBDataHandle,
          typename Data::mapped_type>
  {
  public:
    typedef typename Data::mapped_type DataType;

  public:

    bool contains (int dim, int codim) const
    {
      assert(dim == dimension);
      return codims_.test(codim);
    }

    bool fixedSize (int dim, int codim) const
    {
      assert(dim == dimension);
      return true;
    }

    template<class Entity>
    size_t size (Entity& entity) const
    {
      return 1;
    }

    template<class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& entity) const
    {
      const auto& id = idSet_.id(entity);
      buff.write(data_.at(id));
    }

    template<class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& entity, size_t)
    {
      const auto& id = idSet_.id(entity);
      buff.read(data_[id]);
    }

    LBDataHandle (const IdSet& idSet, Data& data, const Codims& codims)
      : idSet_(idSet)
      , data_(data)
      , codims_(codims)
    {}

  private:
    const IdSet& idSet_;
    Data& data_;
    const Codims codims_;
  };

  template<typename=void>
  static Codims toBitset(Codims codims = {})
    { return codims; }

  template<int codim, int... codimensions, typename=void>
  static Codims toBitset(Codims codims = {})
    { return toBitset<codimensions...>(codims.set(codim)); }

  template<typename=void>
  static void fillVector(const GridView&, Data&) {}

  template<int codim, int... codimensions, typename=void>
  static void fillVector(const GridView& gv, Data& data)
  {
    const auto& idSet = gv.grid().localIdSet();

    for (const auto& entity : entities(gv, Dune::Codim<codim>(),
                                       Dune::Partitions::interiorBorder)) {
      const auto& id = idSet.id(entity);

      // assign the position of the entity to the entry in the vector
      data[id] = entity.geometry().center();
    }

    fillVector<codimensions...>(gv, data);
  }

  template<typename=void>
  static bool checkVector(const GridView&, const Data&)
    { return true; }

  template<int codim, int... codimensions, typename=void>
  static bool checkVector(const GridView& gv, const Data& data)
  {
    const auto& idSet = gv.grid().localIdSet();

    for (const auto& entity : entities(gv, Dune::Codim<codim>(),
                                       Dune::Partitions::interiorBorder)) {
      const auto& id = idSet.id(entity);
      const auto& commPos = data.at(id);
      const auto& realPos = entity.geometry().center();

      // compare the position with the balanced data
      for (int k = 0; k < dimension; k++)
      {
        if (Dune::FloatCmp::ne(commPos[k], realPos[k]))
        {
          DUNE_THROW(Dune::ParallelError,
                     gv.comm().rank() << ": position " << realPos
                                      << " does not coincide with communicated data "
                                      << commPos);
        }
      }
    }

    return checkVector<codimensions...>(gv, data);
  }

public:
  template<int... codimensions>
  static void test(Grid& grid)
  {
    const Codims codims = toBitset<codimensions...>();
    const auto& gv = grid.leafGridView();

    // define the vector containing the data to be balanced
    Data data;

    LBDataHandle dataHandle(grid.localIdSet(), data, codims);

    // fill the data vector
    fillVector<codimensions...>(gv, data);

    // balance the grid and the data
    grid.loadBalance(dataHandle);

    // check for correctness
    checkVector<codimensions...>(gv, data);

    dverb << gv.comm().rank()
              << ": load balancing with data was successful." << std::endl;
  }
};

template<typename Grid>
std::shared_ptr<Grid>
setupGrid(bool simplexGrid, bool localRefinement, int refinementDim, bool refineUpperPart, int numCellsPerDim = 4)
{
  const static int dim = Grid::dimension;
  StructuredGridFactory<Grid> structuredGridFactory;

  Dune::FieldVector<double,dim> lowerLeft(0);
  Dune::FieldVector<double,dim> upperRight(1);
  std::array<unsigned int, dim> numElements;
  std::fill(numElements.begin(), numElements.end(), numCellsPerDim);
  if (simplexGrid)
    return structuredGridFactory.createSimplexGrid(lowerLeft, upperRight, numElements);
  else
    return structuredGridFactory.createCubeGrid(lowerLeft, upperRight, numElements);
}

template<int dim, int... codimensions>
void testLoadBalance(bool simplexGrid, bool localRefinement, int refinementDim, bool refineUpperPart)
{
  using Grid = UGGrid<dim>;
  auto grid = setupGrid<Grid>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
  LoadBalance<Grid>::template test<codimensions...>(*grid);
  // LoadBalance<Grid>::template test<codimensions...>(*grid);
}

template<int dim>
void testParallelUG(bool simplexGrid, bool localRefinement, int refinementDim, bool refineUpperPart)
{
  std::cout << "Testing parallel UGGrid for " << dim << "D\n";

  ////////////////////////////////////////////////////////////
  // Create uniform grid on rank 0
  ////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;
  auto grid = setupGrid<GridType>(simplexGrid, localRefinement, refinementDim, refineUpperPart);

  //////////////////////////////////////////////////////
  // Distribute the grid
  //////////////////////////////////////////////////////
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
  const LeafGV&  leafGridView = grid->leafGridView();
  const LevelGV& level0GridView = grid->levelGridView(0);

  dverb << "LevelGridView for level 0 has " << level0GridView.size(0)
            << " elements "  << level0GridView.size(dim - 1)
            << " edges and " << level0GridView.size(dim)
            << " nodes.\n";
  dverb << "LeafGridView has " << leafGridView.size(0)
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

  // Test communication
  std::map<InterfaceType, std::set<PartitionType> > sendingPartitions;
  sendingPartitions[InteriorBorder_InteriorBorder_Interface] = {InteriorEntity, BorderEntity};
  sendingPartitions[InteriorBorder_All_Interface] = {InteriorEntity, BorderEntity};
  sendingPartitions[Overlap_OverlapFront_Interface] = {OverlapEntity};
  sendingPartitions[Overlap_All_Interface] = {OverlapEntity};
  sendingPartitions[All_All_Interface] = {InteriorEntity, BorderEntity, OverlapEntity, FrontEntity, GhostEntity};

  std::map<InterfaceType, std::set<PartitionType> > receivingPartitions;
  receivingPartitions[InteriorBorder_InteriorBorder_Interface] = {InteriorEntity, BorderEntity};
  receivingPartitions[InteriorBorder_All_Interface] = {InteriorEntity, BorderEntity, OverlapEntity, FrontEntity, GhostEntity};
  receivingPartitions[Overlap_OverlapFront_Interface] = {OverlapEntity, FrontEntity};
  receivingPartitions[Overlap_All_Interface] = {InteriorEntity, BorderEntity, OverlapEntity, FrontEntity, GhostEntity};
  receivingPartitions[All_All_Interface] = {InteriorEntity, BorderEntity, OverlapEntity, FrontEntity, GhostEntity};

  // Test all communication interfaces
  for (auto&& communicationInterface : {InteriorBorder_InteriorBorder_Interface,
                                        InteriorBorder_All_Interface,
                                        Overlap_OverlapFront_Interface,
                                        Overlap_All_Interface,
                                        All_All_Interface})
  {
    // The variable codimSet encodes a set of codimensions as a bitset.
    // We loop over all possible sets.
    for (std::size_t codimSet=0; codimSet<(1<<(dim+1)); codimSet++)
    {
      // TODO: Test BackwardCommunication, too!
      testCommunication<LevelGV>(level0GridView, std::bitset<dim+1>(codimSet),
                                 communicationInterface, ForwardCommunication,
                                 sendingPartitions[communicationInterface],
                                 receivingPartitions[communicationInterface]);

      testCommunication<LeafGV>(leafGridView, std::bitset<dim+1>(codimSet),
                                communicationInterface, ForwardCommunication,
                                sendingPartitions[communicationInterface],
                                receivingPartitions[communicationInterface]);
    }
  }

  ////////////////////////////////////////////////////
  //  Refine globally and test again
  ////////////////////////////////////////////////////

  if (!localRefinement)
  {
    grid->globalRefine(1);
  }
  else
  {
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

  // Test all communication interfaces
  for (auto&& communicationInterface : {InteriorBorder_InteriorBorder_Interface,
                                        InteriorBorder_All_Interface,
                                        Overlap_OverlapFront_Interface,
                                        Overlap_All_Interface,
                                        All_All_Interface})
  {
    // The variable codimSet encodes a set of codimensions as a bitset.
    // We loop over all possible sets.
    for (std::size_t codimSet=0; codimSet<(1<<(dim+1)); codimSet++)
    {
      // TODO: Test BackwardCommunication, too!
      for (int i=0; i<=grid->maxLevel(); i++)
        testCommunication<LevelGV>(grid->levelGridView(i), std::bitset<dim+1>(codimSet),
                                communicationInterface, Dune::ForwardCommunication,
                                sendingPartitions[communicationInterface],
                                receivingPartitions[communicationInterface]);

      testCommunication<LeafGV>(leafGridView, std::bitset<dim+1>(codimSet),
                                communicationInterface, Dune::ForwardCommunication,
                                sendingPartitions[communicationInterface],
                                receivingPartitions[communicationInterface]);
    }
  }
}

template<class GridView>
bool isViewWithinBounds(const GridView& gridView,
                        const Dune::FieldVector<double, GridView::dimensionworld>& lowerLeft,
                        const Dune::FieldVector<double, GridView::dimensionworld>& upperRight)
{
  for (const auto& element : elements(gridView, Dune::Partitions::interior))
  {
    const auto center = element.geometry().center();
    for (int dimIdx = 0; dimIdx < int(GridView::dimensionworld); ++dimIdx)
      if (center[dimIdx] < lowerLeft[dimIdx] || center[dimIdx] > upperRight[dimIdx])
        return false;
  }
  return true;
}

template<class GridView, std::enable_if_t<GridView::dimension==2, int> = 0>
void checkDefaultBisectionPartition(const GridView& gridView, int numProcs, int rank)
{
  // check if all partitions exist (we don't care which processor has which partition)
  int part = -1;
  switch (numProcs){
  case 1:
    if (isViewWithinBounds(gridView, {0,0}, {1,1})) part = 0;
    break;
  case 2:
    if (isViewWithinBounds(gridView, {0,0}, {0.5,1})) part = 0;
    if (isViewWithinBounds(gridView, {0.5,0}, {1,1})) part = 1;
    break;
  case 3:
    if (isViewWithinBounds(gridView, {0,0}, {1.0/3.0,1})) part = 0;
    if (isViewWithinBounds(gridView, {1.0/3.0,0}, {1,0.5})) part = 1;
    if (isViewWithinBounds(gridView, {1.0/3.0,0.5}, {1,1})) part = 2;
    break;
  case 4:
    if (isViewWithinBounds(gridView, {0,0}, {0.5,0.5})) part = 0;
    if (isViewWithinBounds(gridView, {0,0.5}, {0.5,1})) part = 1;
    if (isViewWithinBounds(gridView, {0.5,0}, {1,0.5})) part = 2;
    if (isViewWithinBounds(gridView, {0.5,0.5}, {1,1})) part = 3;
    break;
  case 8:
    if (isViewWithinBounds(gridView, {0,0}, {0.25,0.5})) part = 0;
    if (isViewWithinBounds(gridView, {0.25,0}, {0.5,0.5})) part = 1;
    if (isViewWithinBounds(gridView, {0.5,0}, {0.75,0.5})) part = 2;
    if (isViewWithinBounds(gridView, {0.75,0}, {1,0.5})) part = 3;
    if (isViewWithinBounds(gridView, {0,0.5}, {0.25,1})) part = 4;
    if (isViewWithinBounds(gridView, {0.25,0.5}, {0.5,1})) part = 5;
    if (isViewWithinBounds(gridView, {0.5,0.5}, {0.75,1})) part = 6;
    if (isViewWithinBounds(gridView, {0.75,0.5}, {1,1})) part = 7;
    break;
  default:
    DUNE_THROW(Dune::NotImplemented, "Test for number of processors: " << numProcs);
  }
  const auto sumPart = gridView.comm().sum(part);
  const auto n = numProcs-1;
  if (sumPart != (n*n + n)/2)
    DUNE_THROW(Dune::Exception, "Invalid partition!");
}

template<class GridView, std::enable_if_t<GridView::dimension==3, int> = 0>
void checkDefaultBisectionPartition(const GridView& gridView, int numProcs, int rank)
{
  int part = -1;
  switch (numProcs){
  case 1:
    if (isViewWithinBounds(gridView, {0,0,0}, {1,1,1})) part = 0;
    break;
  case 2:
    if (isViewWithinBounds(gridView, {0,0,0}, {0.5,1,1})) part = 0;
    if (isViewWithinBounds(gridView, {0.5,0,0}, {1,1,1})) part = 1;
    break;
  case 3:
    if (isViewWithinBounds(gridView, {0,0,0}, {1.0/3.0,1,1})) part = 0;
    if (isViewWithinBounds(gridView, {1.0/3.0,0,0}, {1,0.5,1})) part = 1;
    if (isViewWithinBounds(gridView, {1.0/3.0,0.5,0}, {1,1,1})) part = 2;
    break;
  case 4:
    if (isViewWithinBounds(gridView, {0,0,0}, {0.5,0.5,1})) part = 0;
    if (isViewWithinBounds(gridView, {0,0.5,0}, {0.5,1,1})) part = 1;
    if (isViewWithinBounds(gridView, {0.5,0,0}, {1,0.5,1})) part = 2;
    if (isViewWithinBounds(gridView, {0.5,0.5,0}, {1,1,1})) part = 3;
  case 8:
    if (isViewWithinBounds(gridView, {0,0,0}, {0.5,0.5,0.5})) part = 0;
    if (isViewWithinBounds(gridView, {0,0.5,0}, {0.5,1,0.5})) part = 1;
    if (isViewWithinBounds(gridView, {0.5,0,0}, {1,0.5,0.5})) part = 2;
    if (isViewWithinBounds(gridView, {0.5,0.5,0}, {1,1,0.5})) part = 3;
    if (isViewWithinBounds(gridView, {0,0,0.5}, {0.5,0.5,1})) part = 4;
    if (isViewWithinBounds(gridView, {0,0.5,0.5}, {0.5,1,1})) part = 5;
    if (isViewWithinBounds(gridView, {0.5,0,0.5}, {1,0.5,1})) part = 6;
    if (isViewWithinBounds(gridView, {0.5,0.5,0.5}, {1,1,1})) part = 7;
  break;
    default:
    DUNE_THROW(Dune::NotImplemented, "Test for number of processors: " << numProcs);
  return;
  }
  const auto sumPart = gridView.comm().sum(part);
  const auto n = numProcs-1;
  if (sumPart != (n*n + n)/2)
    DUNE_THROW(Dune::Exception, "Invalid partition!");
}

template<int dim>
void testDefaultLoadBalanceStructuredCube(const Dune::MPIHelper& mpiHelper)
{
  if (mpiHelper.rank() == 0)
    std::cout << "Testing default load balancer for structured cube grid in " << dim << "D\n";

  auto grid = setupGrid<UGGrid<dim>>(false, false, 0, false, mpiHelper.size()*2);
  grid->loadBalance();
  grid->globalRefine(1);
  const auto& gridView = grid->leafGridView();
  checkDefaultBisectionPartition(gridView, mpiHelper.size(), mpiHelper.rank());
}

int main (int argc , char **argv) try
{
  // initialize MPI, finalize is done automatically on exit
  auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

#if !defined(ModelP)
  if (mpiHelper.size() > 1) {
    if (mpiHelper.rank() == 0)
      std::cout << "SKIPPED: test-parallel-ug requires a parallel version of dune-uggrid\n";
    return 77;
  }
#endif

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
        for (const int refinementDim : {0,1}) {
          testParallelUG<2>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<2>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<2, 0>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<2, 2>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<2, 0, 2>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
        }
        for (const int refinementDim : {0,1,2}) {
          testParallelUG<3>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<3>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<3, 0>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<3, 3>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
          testLoadBalance<3, 0, 3>(simplexGrid, localRefinement, refinementDim, refineUpperPart);
        }
      }
    }
  }

  // test default bisection load balancer for case with known partition
  testDefaultLoadBalanceStructuredCube<2>(mpiHelper);
  testDefaultLoadBalanceStructuredCube<3>(mpiHelper);

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
