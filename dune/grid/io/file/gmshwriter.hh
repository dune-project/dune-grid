// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_GMSHWRITER_HH
#define DUNE_GRID_IO_FILE_GMSHWRITER_HH

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

  /**
     \ingroup Gmsh

     \brief Write Gmsh mesh file

     Write a grid using the given GridView as an ASCII Gmsh file of version 2.0.

     If the grid contains an element type not supported by gmsh an IOError exception is thrown.

     All grids in a gmsh file live in three-dimensional Euclidean space. If the world dimension
     of the grid type that you are writing is less than three, the remaining coordinates are
     set to zero.
   */
  template <class GridView>
  class GmshWriter
  {
  private:
    const GridView gv;
    int precision;

    static const unsigned int dim = GridView::dimension;
    static const unsigned int dimWorld = GridView::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshWriter requires dimWorld <= 3." );

    /** \brief Returns index of i-th vertex of an element, plus 1 (for gmsh numbering) */
    template<typename Entity>
    std::size_t nodeIndexFromEntity(const Entity& entity, int i) const {
      return gv.indexSet().subIndex(entity, i, dim)+1;
    }

    /** \brief Translate GeometryType to corresponding Gmsh element type number
      * \throws IOError if there is no equivalent type in Gmsh
      */
    static std::size_t translateDuneToGmshType(const GeometryType& type) {
      std::size_t element_type;

      if (type.isLine())
        element_type = 1;
      else if (type.isTriangle())
        element_type = 2;
      else if (type.isQuadrilateral())
        element_type = 3;
      else if (type.isTetrahedron())
        element_type = 4;
      else if (type.isHexahedron())
        element_type = 5;
      else if (type.isPrism())
        element_type = 6;
      else if (type.isPyramid())
        element_type = 7;
      else if (type.isVertex())
        element_type = 15;
      else
        DUNE_THROW(Dune::IOError, "GeometryType " << type << " is not supported by gmsh.");

      return element_type;
    }

    /** \brief Writes all the elements of a grid line by line
     *
     * Each line has the format
     *    element-number element-type number-of-tags <tags> node-number-list
     * Counting of the element numbers starts by "1".
     *
     * If `physicalEntities` is not empty, each element has a tag representing its physical id.
     *
     * If `physicalBoundaries` is not empty, also the boundaries are written to the file with
     * the corresponding physical value.
     *
     * The physicalBoundaries vector need to be sorted according to the interesection
     * boundary segment index.
     */
    void outputElements(std::ofstream& file, const std::vector<int>& physicalEntities, const std::vector<int>& physicalBoundaries) const {
      MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(gv, mcmgElementLayout());
      std::size_t counter(1);
      for (const auto& entity : elements(gv)) {
        // Check whether the type is compatible. If not, close file and rethrow exception.
        try {
          std::size_t element_type = translateDuneToGmshType(entity.type());

          file << counter << " " << element_type;
          // If present, set the first tag to the physical entity
          if (!physicalEntities.empty())
            file << " " << 1 << " " << physicalEntities[elementMapper.index(entity)];
          else
            file << " " << 0; // "0" for "I do not use any tags."

          // Output list of nodes.
          // 3, 5 and 7 got different vertex numbering compared to Dune
          if (3 == element_type)
            file << " "
                 << nodeIndexFromEntity(entity, 0) << " " << nodeIndexFromEntity(entity, 1) << " "
                 << nodeIndexFromEntity(entity, 3) << " " << nodeIndexFromEntity(entity, 2);
          else if (5 == element_type)
            file << " "
                 << nodeIndexFromEntity(entity, 0) << " " << nodeIndexFromEntity(entity, 1) << " "
                 << nodeIndexFromEntity(entity, 3) << " " << nodeIndexFromEntity(entity, 2) << " "
                 << nodeIndexFromEntity(entity, 4) << " " << nodeIndexFromEntity(entity, 5) << " "
                 << nodeIndexFromEntity(entity, 7) << " " << nodeIndexFromEntity(entity, 6);
          else if (7 == element_type)
            file << " "
                 << nodeIndexFromEntity(entity, 0) << " " << nodeIndexFromEntity(entity, 1) << " "
                 << nodeIndexFromEntity(entity, 3) << " " << nodeIndexFromEntity(entity, 2) << " "
                 << nodeIndexFromEntity(entity, 4);
          else {
            for (int k = 0; k < entity.geometry().corners(); ++k)
              file << " " << nodeIndexFromEntity(entity, k);
          }
          ++counter;

          file << std::endl;

          // Write boundaries
          if (!physicalBoundaries.empty()) {
            auto refElement = referenceElement<typename GridView::ctype,dim>(entity.type());
            for(const auto& intersection : intersections(gv, entity)) {
              if(intersection.boundary()) {
                const auto faceLocalIndex(intersection.indexInInside());
                file << counter << " " << translateDuneToGmshType(intersection.type())
                  << " " << 1 << " " << physicalBoundaries[intersection.boundarySegmentIndex()];
                for (int k = 0; k < intersection.geometry().corners(); ++k)
                {
                  const auto vtxLocalIndex(refElement.subEntity(faceLocalIndex, 1, k, dim));
                  file << " " << nodeIndexFromEntity(entity, vtxLocalIndex);
                }
                ++counter;
                file << std::endl;
              }
            }
          }

        } catch(Exception& e) {
          file.close();
          throw;
        }
      }
    }


    /** \brief Writes all the vertices of a grid line by line
     *
     * Each line has the format
     *  node-number x-coord y-coord z-coord
     * The node-numbers will most certainly not have the arrangement "1, 2, 3, ...".
     */
    void outputNodes(std::ofstream& file) const {
      for (const auto& vertex : vertices(gv)) {
        const auto globalCoord = vertex.geometry().center();
        const auto nodeIndex = gv.indexSet().index(vertex)+1; // Start counting indices by "1".

        if (1 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << 0 << " " << 0 << std::endl;
        else if (2 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << globalCoord[1] << " " << 0 << std::endl;
        else // (3 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << globalCoord[1] << " " << globalCoord[2] << std::endl;
      }
    }

  public:
    /**
     * \brief Constructor expecting GridView of Grid to be written.
     * \param gridView GridView that will be written.
     * \param numDigits Number of digits to use.
     */
    GmshWriter(const GridView& gridView, int numDigits=6) : gv(gridView), precision(numDigits) {}

    /**
     * \brief Set the number of digits to be used when writing the vertices. By default is 6.
     * \brief numDigits Number of digits to use.
     */
    void setPrecision(int numDigits) {
      precision = numDigits;
    }

    /**
     * \brief Write given grid in Gmsh 2.0 compatible ASCII file.
     * \param fileName Path of file. This method does not attach a ".msh"-extension by itself.
     * \param physicalEntities Physical entities for each element (optional).
     * \param physicalBoundaries Physical boundaries (optional).
     *
     * Opens the file with given name and path, stores the element data of the grid
     * and closes the file when done.
     *
     * If the optional parameter `physicalEntities` is provided, each element is written with
     * a tag representing its physical id.
     *
     * If the optional parameter `physicalBoundaries` is provided, also the boundaries
     * are written on file with the corresponding physical value.
     *
     * The physicalBoundaries vector need to be sorted according to the interesection boundary
     * segment index.
     *
     * Throws an IOError if file could not be opened or an unsupported element type is
     * encountered.
     */
    void write(const std::string& fileName,
               const std::vector<int>& physicalEntities=std::vector<int>(),
               const std::vector<int>& physicalBoundaries=std::vector<int>()) const {
      // Open file
      std::ofstream file(fileName.c_str());
      if (!file.is_open())
        DUNE_THROW(Dune::IOError, "Could not open " << fileName << " with write access.");

      // Set precision
      file << std::setprecision( precision );

      // Output Header
      file << "$MeshFormat" << std::endl
           << "2.0 0 " << sizeof(double) << std::endl // "2.0" for "version 2.0", "0" for ASCII
           << "$EndMeshFormat" << std::endl;

      // Output Nodes
      file << "$Nodes" << std::endl
           << gv.size(dim) << std::endl;

      outputNodes(file);

      file << "$EndNodes" << std::endl;

      // Output Elements;
      int boundariesSize(0);
      if(!physicalBoundaries.empty())
        for(const auto& entity : elements(gv))
          for(const auto& intersection : intersections(gv, entity))
            if(intersection.boundary())
              ++boundariesSize;

      file << "$Elements" << std::endl
           << gv.size(0) + boundariesSize<< std::endl;

      outputElements(file, physicalEntities, physicalBoundaries);

      file << "$EndElements" << std::endl;
    }

  };

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_GMSHWRITER_HH
