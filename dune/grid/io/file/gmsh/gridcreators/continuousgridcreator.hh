// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_GRIDCREATORS_CONTINUOUSGRIDCREATOR_HH
#define DUNE_GRID_IO_FILE_GMSH_GRIDCREATORS_CONTINUOUSGRIDCREATOR_HH

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/io/file/gmsh/types.hh>
#include <dune/grid/io/file/gmsh/gridcreatorinterface.hh>
#include <dune/grid/io/file/gmsh/utility/errors.hh>

namespace Dune::Impl::Gmsh
{
  // Create a grid where the input points and connectivity is already
  // connected correctly.
  template <class Grid>
  struct ContinuousGridCreator
    : public GridCreatorInterface<Grid, ContinuousGridCreator<Grid> >
  {
    using Super = GridCreatorInterface<Grid, ContinuousGridCreator<Grid> >;
    using GlobalCoordinate = typename Super::GlobalCoordinate;
    using Nodes = std::vector<GlobalCoordinate>;


  public:
    using Super::Super;
    using Super::factory;

    template <class NodeAttributes>
    void insertVerticesImpl (std::size_t numNodes,
                             std::pair<std::size_t,std::size_t> nodeTagRange,
                             std::vector<NodeAttributes> const& entityBlocks)
    {
      vertexMap_.resize(nodeTagRange.second - nodeTagRange.first + 1);
      vertexShift_ = nodeTagRange.first;
      nodes_.resize(numNodes);
      GlobalCoordinate p;
      size_t vertexIndex = 0;

      for (auto const& entityBlock : entityBlocks) {
        for (auto const& node : entityBlock.nodes) {
          for (std::size_t j = 0; j < p.size(); ++j)
            p[j] = node.xyz[j];
          nodes_[vertexIndex] = p;
          vertexMap_[node.tag - vertexShift_] = vertexIndex++;
        }
      }
    }

    template <class ElementAttributes, class BoundaryEntities>
    void insertElementsImpl (std::size_t /*numElements*/,
                             std::pair<std::size_t,std::size_t> /*elementTagRange*/,
                             std::vector<ElementAttributes> const& entityBlocks,
                             BoundaryEntities const& /*boundaryEntities*/)
    {
      std::vector<unsigned int> connectivity;
      std::size_t cornerIndex = 0;
      std::vector<std::int64_t> cornerVertices(nodes_.size(), -1);

      for (auto const& entityBlock : entityBlocks) {
        if (entityBlock.entityDim < Grid::dimension-1)
          continue;

        auto type = gmshNumberToGeometryType(entityBlock.elementType);
        CellType cell{type};

        if (entityBlock.entityDim == Grid::dimension) {   //element
          auto refElem = referenceElement<double,Grid::dimension>(cell.type());
          connectivity.resize(refElem.size(Grid::dimension));

          for (auto const& element : entityBlock.elements) {
            GMSH4_ASSERT(element.nodes.size() >= connectivity.size());
            for (std::size_t j = 0; j < connectivity.size(); ++j) {
              auto index = vertexMap_[element.nodes[j] - vertexShift_];
              auto& vertex = cornerVertices.at(index);
              if (vertex < 0) {
                factory().insertVertex(nodes_.at(index));
                vertex = cornerIndex++;
              }
              connectivity[cell.gmshVertexToDuneVertex(j)] = vertex;
            }

            factory().insertElement(cell.type(), connectivity);
          }
        }
      }
      nodes_.clear();
    }

  private:
    /// All point coordinates including the higher-order Lagrange points
    Nodes nodes_;

    std::vector<std::size_t> vertexMap_;
    std::size_t vertexShift_ = 0;
  };

  // deduction guides
  template <class Grid>
  ContinuousGridCreator(GridFactory<Grid>&)
  -> ContinuousGridCreator<Grid>;

} // end namespace Dune::Impl::Gmsh

#endif
