// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <memory>

#include <dune/common/stdstreams.hh>

#include <dune/grid/uggrid/uggridfactory.hh>
#include "boundaryextractor.hh"

#include <dune-uggrid-config.hh>
#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
#  include <dune/uggrid/parallel/ppif/ppifcontext.hh>
#endif

namespace Dune {

/* The following three methods are the ones that UG calls to know about the geometry
   of the boundary.  UG expects one static method for each coarse grid boundary segment.
   In order to let separate BoundarySegment objects handle the calls we let UG pass
   them along as 'data'.

   \param data The BoundarySegment object
   \param param Local coordinates on the boundary segment
   \param result The world coordinates of the corresponding boundary point
 */

static int boundarySegmentWrapper2d(void *data, double *param, FieldVector<double,2>& result)
{
  const BoundarySegment<2>* boundarySegment = static_cast<BoundarySegment<2>*>(data);

  result = (*boundarySegment)(*((FieldVector<double,1>*)param));

  return 0;
}

static int boundarySegmentWrapper3dTriangle(void *data, double *param, FieldVector<double,3>& result)
{
  const BoundarySegment<3>* boundarySegment = static_cast<BoundarySegment<3>*>(data);

  // UG produces coordinates in the triangle (0,0) -- (1,0) -- (1,1),
  // but we want them to be in (0,0) -- (1,0) -- (0,1)
  FieldVector<double, 2> local;
  local[0] = param[0]-param[1];
  local[1] = param[1];

  result = (*boundarySegment)(local);

  return 0;
}

static int boundarySegmentWrapper3dQuad(void *data, double *param, FieldVector<double,3>& result)
{
  const BoundarySegment<3>* boundarySegment = static_cast<BoundarySegment<3>*>(data);

  result = (*boundarySegment)(*((FieldVector<double,2>*)param));

  return 0;
}



template <int dimworld>
GridFactory<UGGrid<dimworld> >::
GridFactory()
{
  grid_ = new UGGrid<dimworld>;

  factoryOwnsGrid_ = true;

  createBegin();
}

template <int dimworld>
GridFactory<UGGrid<dimworld> >::
GridFactory(UGGrid<dimworld>* grid)
{
  grid_ = grid;

  factoryOwnsGrid_ = false;

  createBegin();
}

template <int dimworld>
GridFactory<UGGrid<dimworld> >::
~GridFactory()
{
  if (grid_ && factoryOwnsGrid_)
    delete grid_;
}

template <int dimworld>
void GridFactory<UGGrid<dimworld> >::
insertVertex(const FieldVector<typename GridFactory<UGGrid<dimworld> >::ctype,dimworld>& pos)
{
  vertexPositions_.push_back(pos);
}

template <int dimworld>
void GridFactory<UGGrid<dimworld> >::
insertElement(const GeometryType& type,
              const std::vector<unsigned int>& vertices)
{
  if (dimworld!=type.dim())
    DUNE_THROW(GridError, "You cannot insert a " << type
                                                 << " into a UGGrid<" << dimworld << ">!");

  int newIdx = elementVertices_.size();

  elementTypes_.push_back(vertices.size());
  for (unsigned int i=0; i<vertices.size(); i++)
    elementVertices_.push_back(vertices[i]);

  if (type.isTriangle()) {
    // Everything alright
    if (vertices.size() != 3)
      DUNE_THROW(GridError, "You have requested to enter a triangle, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isQuadrilateral()) {

    if (vertices.size() != 4)
      DUNE_THROW(GridError, "You have requested to enter a quadrilateral, but you"
                 << " have provided " << vertices.size() << " vertices!");

    // DUNE and UG numberings differ --> reorder the vertices
    elementVertices_[newIdx+2] = vertices[3];
    elementVertices_[newIdx+3] = vertices[2];

  } else if (type.isTetrahedron()) {

    if (vertices.size() != 4)
      DUNE_THROW(GridError, "You have requested to enter a tetrahedron, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isPyramid()) {

    if (vertices.size() != 5)
      DUNE_THROW(GridError, "You have requested to enter a pyramid, but you"
                 << " have provided " << vertices.size() << " vertices!");

    // DUNE and UG numberings differ --> reorder the vertices
    elementVertices_[newIdx+2] = vertices[3];
    elementVertices_[newIdx+3] = vertices[2];

  } else if (type.isPrism()) {

    if (vertices.size() != 6)
      DUNE_THROW(GridError, "You have requested to enter a prism, but you"
                 << " have provided " << vertices.size() << " vertices!");

  } else if (type.isHexahedron()) {

    if (vertices.size() != 8)
      DUNE_THROW(GridError, "You have requested to enter a hexahedron, but you"
                 << " have provided " << vertices.size() << " vertices!");

    // DUNE and UG numberings differ --> reorder the vertices
    elementVertices_[newIdx+2] = vertices[3];
    elementVertices_[newIdx+3] = vertices[2];
    elementVertices_[newIdx+6] = vertices[7];
    elementVertices_[newIdx+7] = vertices[6];

  } else {
    DUNE_THROW(GridError, "You cannot insert a " << type
                                                 << " into a UGGrid<" << dimworld << ">!");
  }

}

template <int dimworld>
void GridFactory<UGGrid<dimworld> >::
insertBoundarySegment(const std::vector<unsigned int>& vertices)
{
  insertBoundarySegment(vertices, std::shared_ptr<BoundarySegment<dimworld> >());
}

template <int dimworld>
void GridFactory<UGGrid<dimworld> >::
insertBoundarySegment(const std::vector<unsigned int>& vertices,
                      const std::shared_ptr<BoundarySegment<dimworld> > &boundarySegment)
{
  std::array<int, dimworld*2-2> segmentVertices;

  for (size_t i=0; i<vertices.size(); i++)
    segmentVertices[i] = vertices[i];

  for (size_t i=vertices.size(); i<dimworld*2-2; i++)
    segmentVertices[i] = -1;

  // DUNE --> UG vertex renumbering for quadrilateral boundary segments (3D only)
  if constexpr (dimworld == 3) {
    if (vertices.size()==4) {
      segmentVertices[2] = vertices[3];
      segmentVertices[3] = vertices[2];
    }
  }

  boundarySegmentVertices_.push_back(segmentVertices);

  // Append boundary segment class to the boundary segment class list,
  // to make sure they aren't deleted before the grid object
  grid_->boundarySegments_.push_back(boundarySegment);
}

template <int dimworld>
std::unique_ptr<UGGrid<dimworld>> GridFactory<UGGrid<dimworld> >::
createGrid()
{
  // Prevent a crash when this method is called twice in a row
  // You never know who may do this...
  if (grid_==nullptr)
    return nullptr;

  // ///////////////////////////////////////////////////////////////////////////////
  //  Communicate the grid information from the master to all other processes
  //
  // Rationale: Parallel Dune grids that are set up by providing a coarse grid
  // on rank 0 and then loadbalancing should not require that the grid be
  // provided to all processes directly. Usually UGGrid does that, because
  // it needs the grid only on rank 0, but the grid _boundary_ on all processes.
  // This patch makes the rank 0 process broadcast the entire coarse
  // grid information from rank 0 to all other processes before
  // starting the actual grid creation.  That way, the user has to
  // supply the coarse grid only on rank 0 (as the specification says),
  // but still the UGGridFactory has all relevant information available
  // on all processes.
  // ///////////////////////////////////////////////////////////////////////////////

  const auto& comm = grid_->comm();

  // Broadcast the vertex positions
  int numVertices = vertexPositions_.size();
  comm.broadcast(&numVertices, 1, 0);
  vertexPositions_.resize(numVertices);
  comm.broadcast(&vertexPositions_[0], vertexPositions_.size(), 0);

  // Broadcast the element corners
  int numElementVertices = elementVertices_.size();
  comm.broadcast(&numElementVertices, 1, 0);
  elementVertices_.resize(numElementVertices);
  comm.broadcast(&elementVertices_[0], elementVertices_.size(), 0);

  // Broadcast the geometry types
  int numElementTypes = elementTypes_.size();
  comm.broadcast(&numElementTypes, 1, 0);
  elementTypes_.resize(numElementTypes);
  comm.broadcast(&elementTypes_[0], elementTypes_.size(), 0);

  // ///////////////////////////////////////////
  //   Extract grid boundary segments
  // ///////////////////////////////////////////
  std::set<UGGridBoundarySegment<dimworld> > boundarySegments;
  typedef typename std::set<UGGridBoundarySegment<dimworld> >::iterator SetIterator;

  BoundaryExtractor::detectBoundarySegments(elementTypes_, elementVertices_, boundarySegments);
  if (boundarySegments.empty())
    DUNE_THROW(GridError, "Couldn't extract grid boundary.");

  std::vector<int> isBoundaryNode;
  BoundaryExtractor::detectBoundaryNodes(boundarySegments, vertexPositions_.size(), isBoundaryNode);

  dverb << boundarySegments.size() << " boundary segments were found!" << std::endl;

  if (grid_->boundarySegments_.size() > boundarySegments.size())
    DUNE_THROW(GridError, "You have supplied " << grid_->boundarySegments_.size()
                                               << " parametrized boundary segments, but the coarse grid has only "
                                               << boundarySegments.size() << " boundary faces!");

  // Count number of nodes on the boundary
  int noOfBNodes = 0;
  for (unsigned int i=0; i<isBoundaryNode.size(); i++) {
    if (isBoundaryNode[i] != -1)
      noOfBNodes++;
  }

  // ///////////////////////////////////////////////////////
  //   UG needs all boundary vertices first.  We set up
  //   up an array to keep track of the reordering.
  // ///////////////////////////////////////////////////////
  std::vector<unsigned int> nodePermutation;

  for (unsigned int i=0; i<isBoundaryNode.size(); ++i) {
    if (isBoundaryNode[i]!=-1)
      nodePermutation.push_back(i);
  }

  for (unsigned int i=0; i<isBoundaryNode.size(); ++i) {
    if (isBoundaryNode[i]==-1)
      nodePermutation.push_back(i);
  }

  // ///////////////////////////////////////////
  //   Create the domain data structure
  // ///////////////////////////////////////////
  grid_->numBoundarySegments_ = boundarySegments.size();

  auto ugDomain = std::make_unique<typename UG_NS<dimworld>::domain>();
  ugDomain->numOfSegments = grid_->numBoundarySegments_;
  ugDomain->numOfCorners = noOfBNodes;

  // ///////////////////////////////////////////
  //   Insert the boundary segments
  // ///////////////////////////////////////////
  unsigned int i;
  for (i=0; i<grid_->boundarySegments_.size(); i++) {

    // Copy the vertices into a C-style array
    // We copy four vertices here even if the segment is a triangle -- it doesn't matter
    int vertices_c_style[dimworld*2-2];
    for (int j=0; j<dimworld*2-2; j++) {
      // For parameterized boundary segments, -1 is used as a sentinel value
      // for unused vertices and must be preserved.
      const auto idx = boundarySegmentVertices_[i][j];
      vertices_c_style[j] = idx == -1 ? -1 : isBoundaryNode[idx];
    }

    if (grid_->boundarySegments_[i]) {


      typedef int (*BndSegFuncPtr)(void *, double *, FieldVector<double,dimworld>&);
      BndSegFuncPtr boundarySegmentWrapper;
      if constexpr (dimworld==2)
        boundarySegmentWrapper = boundarySegmentWrapper2d;
      else {
        if (vertices_c_style[3]==-1)      // triangular face
          boundarySegmentWrapper = boundarySegmentWrapper3dTriangle;
        else
          boundarySegmentWrapper = boundarySegmentWrapper3dQuad;
      }

      // Actually create the segment
      ugDomain->boundarySegments.emplace_back(i,  // Index of the segment
                                              vertices_c_style,
                                              boundarySegmentWrapper,
                                              grid_->boundarySegments_[i].get());

    } else {       // No explicit boundary parametrization has been given. We insert a linear segment

      int numVertices = (dimworld==2)
                        ? 2
                        : ((boundarySegmentVertices_[i][3]==-1) ? 3 : 4);

      std::array<FieldVector<double,dimworld>, 2*(dimworld-1)> segmentCoordinates;
      for (int j=0; j<numVertices; j++)
        segmentCoordinates[j] = vertexPositions_[boundarySegmentVertices_[i][j]];

      ugDomain->linearSegments.emplace_back(i,                  // id of segment
                                            numVertices,        // Number of corners
                                            vertices_c_style,
                                            segmentCoordinates);
    }

    // /////////////////////////////////////////////////////////////////////
    //   Remove this segment from the set of computed boundary segments.
    // /////////////////////////////////////////////////////////////////////

    UGGridBoundarySegment<dimworld> thisSegment;
    /** \todo Not nice: we need to copy because the array types are different */
    for (int j=0; j<2*dimworld-2; j++)
      thisSegment[j] = boundarySegmentVertices_[i][j];

    SetIterator boundaryElementFace = boundarySegments.find(thisSegment);

    if (boundaryElementFace==boundarySegments.end())
      DUNE_THROW(GridError, "You have provided a boundary parametrization for"
                 << " a segment which is not boundary segment in the grid!");

    // Everything is fine.  Delete the element face from the list to mark that it has been properly handled
    boundarySegments.erase(thisSegment);
  }


  // ///////////////////////////////////////////////////////////////////////
  //   The boundary segments remaining in the std::set boundarySegments
  //   have not been provided with an explicit parametrization.  They are
  //   inserted into the domain as straight boundary segments.
  // ///////////////////////////////////////////////////////////////////////
  SetIterator it = boundarySegments.begin();

  for (; it != boundarySegments.end(); ++it, ++i) {

    const UGGridBoundarySegment<dimworld>& thisSegment = *it;

    // Copy the vertices into a C-style array
    int vertices_c_style[4];

    for (int j=0; j<thisSegment.numVertices(); j++)
      vertices_c_style[j] = isBoundaryNode[thisSegment[j]];

    std::array<FieldVector<double,dimworld>, 2*(dimworld-1)> segmentCoordinates;
    for (int j=0; j<thisSegment.numVertices(); j++)
      segmentCoordinates[j] = vertexPositions_[thisSegment[j]];

    ugDomain->linearSegments.emplace_back(i,                  // id of segment
                                          thisSegment.numVertices(),        // Number of corners
                                          vertices_c_style,
                                          segmentCoordinates);
  }


  /////////////////////////////////////////////////
  //   Set up a UG multigrid data structure
  /////////////////////////////////////////////////

  std::string MultiGridName = grid_->name_;
  std::string FormatName = "DuneFormat" + std::to_string(dimworld) + "d";

  grid_->bvp_->Domain = std::move(ugDomain);

  // Make sure there is no old multigrid object with the same name.
  // TODO: Can this happen at all?
  typename UG_NS<dimworld>::MultiGrid *theMG = UG_NS<dimworld>::GetMultigrid(MultiGridName.c_str());
  if (theMG && DisposeMultiGrid(theMG)!=0)
      DUNE_THROW(GridError, "Closing old multigrid failed!");

  // Create the actual multigrid data structure
  std::shared_ptr<PPIF::PPIFContext> ppifContext;
#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
  ppifContext = std::make_shared<PPIF::PPIFContext>(grid_->comm());
#endif
  grid_->multigrid_ = UG_NS<dimworld>::CreateMultiGrid(MultiGridName.c_str(),grid_->bvp_,FormatName.c_str(),true,true, ppifContext);
  if (grid_->multigrid_==nullptr)
    DUNE_THROW(GridError, "Could not create multigrid");

  // ///////////////////////////////////////////////////////////////
  // If we are in a parallel setting and we are _not_ the master
  // process we can stop here.
  // ///////////////////////////////////////////////////////////////
#if ModelP and DUNE_UGGRID_HAVE_PPIFCONTEXT
  const bool isMaster = grid_->multigrid_->ppifContext().isMaster();
#else
  const bool isMaster = (PPIF::me == PPIF::master);
#endif
  if (not isMaster) {
    // Complete the UG-internal grid data structure even if we are
    // not the master process. (CreateAlgebra communicates via MPI
    // so we would be out of sync if we don't do this here...)
    if (CreateAlgebra(grid_->multigrid_) != UG_NS<dimworld>::GM_OK)
      DUNE_THROW(IOError, "Call of 'UG::D" << dimworld << "::CreateAlgebra' failed!");

    // Create an empty levelIndexSet for level 0.
    // Even though it's empty until the first load balancing
    // it should still be there in case someone wants to access it.
    grid_->levelIndexSets_.resize(1);
    grid_->levelIndexSets_[0] = std::make_shared<UGGridLevelIndexSet<const UGGrid<dimworld> > >();

    /* here all temp memory since CreateMultiGrid is released */
    ReleaseTmpMem(grid_->multigrid_->theHeap, grid_->multigrid_->MarkKey);
    grid_->multigrid_->MarkKey = 0;

    // ///////////////////////////////////////////////////
    // hand over the grid and delete the member pointer
    // ///////////////////////////////////////////////////

    UGGrid<dimworld>* tmp = grid_;
    grid_ = nullptr;
    return std::unique_ptr<UGGrid<dimworld>>(tmp);
  }

  // ////////////////////////////////////////////////
  //   Actually insert the interior vertices
  // ////////////////////////////////////////////////
  int nodeCounter = noOfBNodes;
  for (size_t i=0; i<vertexPositions_.size(); i++) {
    if (isBoundaryNode[i] != -1)
      continue;
    if (UG_NS<dimworld>::InsertInnerNode(grid_->multigrid_->grids[0], &((vertexPositions_[i])[0])) == nullptr)
      DUNE_THROW(GridError, "Inserting a vertex into UGGrid failed!");

    isBoundaryNode[i] = nodeCounter++;
  }

  vertexPositions_.resize(0);

  // ////////////////////////////////////////////////
  //   Actually insert all the elements
  // ////////////////////////////////////////////////

  std::vector<typename UG_NS<dimworld>::Node*> nodePointers(isBoundaryNode.size());
  for (typename UG_NS<dimworld>::Node* theNode=UG_NS<dimworld>::FirstNode(grid_->multigrid_->grids[0]); theNode!=nullptr; theNode=theNode->succ)
    nodePointers[theNode->id] = theNode;

  int idx = 0;
  for (size_t i=0; i<elementTypes_.size(); i++) {

    std::vector<typename UG_NS<dimworld>::Node*> vertices(elementTypes_[i]);
    for (size_t j=0; j<elementTypes_[i]; j++)
      vertices[j] = nodePointers[isBoundaryNode[elementVertices_[idx++]]];

    if (InsertElement(grid_->multigrid_->grids[0], elementTypes_[i],&(vertices[0]),nullptr,nullptr,nullptr)==nullptr)
      DUNE_THROW(GridError, "Inserting element into UGGrid failed!");
  }

  // Not needed any more
  elementTypes_.resize(0);
  elementVertices_.resize(0);

  // Complete the UG-internal grid data structure
  if (CreateAlgebra(grid_->multigrid_) != UG_NS<dimworld>::GM_OK)
    DUNE_THROW(IOError, "Call of 'UG::D" << dimworld << "::CreateAlgebra' failed!");

  /* here all temp memory since CreateMultiGrid is released */
  ReleaseTmpMem(grid_->multigrid_->theHeap, grid_->multigrid_->MarkKey);
  grid_->multigrid_->MarkKey = 0;

  // Set the local indices
  grid_->setIndices(true, &nodePermutation);

  // Clear refinement flags
  for (const auto& element : elements(grid_->levelGridView(0)))
    UG_NS<dimworld>::WriteCW(element.impl().target_, UG_NS<dimworld>::NEWEL_CE, 0);


  // ///////////////////////////////////////////////////
  // hand over the grid and delete the member pointer
  // ///////////////////////////////////////////////////

  UGGrid<dimworld>* tmp = grid_;
  grid_ = nullptr;
  return std::unique_ptr<UGGrid<dimworld>>(tmp);
}

template <int dimworld>
void GridFactory<UGGrid<dimworld> >::
createBegin()
{
  // ///////////////////////////////////////////////////////
  //   Clean up existing grid structure if there is one
  // ///////////////////////////////////////////////////////

  // Delete the UG multigrid if there is one (== createEnd() has already
  // been called once for this object)
  if (grid_->multigrid_) {
    // Set UG's currBVP variable to the BVP corresponding to this
    // grid.  This is necessary if we have more than one UGGrid in use.
    // DisposeMultiGrid will crash if we don't do this
    //UG_NS<dim>::Set_Current_BVP(grid_->multigrid_->theBVP);
    // set the multigrid's bvp pointer to NULL to make sure the BVP
    // is not deleted
    grid_->multigrid_->theBVP = nullptr;
    UG_NS<dimworld>::DisposeMultiGrid(grid_->multigrid_);
    grid_->multigrid_ = nullptr;
  }

  // Delete levelIndexSets if there are any
  grid_->levelIndexSets_.resize(0);

  // //////////////////////////////////////////////////////////
  //   Clear all buffers used during coarse grid creation
  // //////////////////////////////////////////////////////////
  grid_->boundarySegments_.resize(0);
  boundarySegmentVertices_.resize(0);
  elementTypes_.resize(0);
  elementVertices_.resize(0);
  vertexPositions_.resize(0);
}




// Explicit template instantiation
template class GridFactory<UGGrid<2> >;
template class GridFactory<UGGrid<3> >;

} /* namespace Dune */
