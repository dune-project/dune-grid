// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_GMSHWRITER_HH
#define DUNE_GRID_IO_FILE_GMSHWRITER_HH


#include <fstream>
#include <iostream>

#include <string>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/grid.hh>


namespace Dune {

  /**
     \ingroup Gmsh

     \brief Write Gmsh mesh file

     Write a grid using the given GridView as an ASCII Gmsh file of version 2.0.

     If the grid contains an element type not supported by gmsh an IOError exception is thrown.

     Boundary segments are ignored.

     All grids in a gmsh file live in three-dimensional Euclidean space. If the world dimension
     of the grid type that you are writing is less than three, the remaining coordinates are
     set to zero.
   */
  template <class GridView>
  class GmshWriter
  {
  private:
    const GridView gv;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    dune_static_assert( (dimWorld <= 3), "GmshWriter requires dimWorld <= 3." );

    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;


    /** \brief Returns index of i-th vertex of an element, plus 1 (for gmsh numbering) */
    size_t nodeIndexFromIterator(const ElementIterator& eIt, int i) const {
      return gv.indexSet().subIndex(*eIt, i, dim)+1;
    }

    /** \brief Translate GeometryType to corresponding Gmsh element type number
      * \throws IOError if there is no equivalent type in Gmsh
      */
    static size_t translateDuneToGmshType(const GeometryType& type) {
      // Probably the non-clever, but hopefully readable way of translating the GeometryTypes to the gmsh types
      size_t element_type;

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
     * Tags are ignored, i.e. number-of-tags is always zero and no tags are printed.
     * node-number-list depends on the type of the given element.
     */
    void outputElements(std::ofstream& file) const {
      ElementIterator eIt    = gv.template begin<0>();
      ElementIterator eEndIt = gv.template end<0>();

      for (size_t i = 1; eIt != eEndIt; ++eIt, ++i) {
        // Check whether the type is compatible. If not, close file and rethrow exception.
        try {
          size_t element_type = translateDuneToGmshType(eIt->type());

          file << i << " " << element_type << " " << 0; // "0" for "I do not use any tags."

          // Output list of nodes.
          // 3, 5 and 7 got different vertex numbering compared to Dune
          if (3 == element_type)
            file << " "
                 << nodeIndexFromIterator(eIt, 0) << " " << nodeIndexFromIterator(eIt, 1) << " "
                 << nodeIndexFromIterator(eIt, 3) << " " << nodeIndexFromIterator(eIt, 2);
          else if (5 == element_type)
            file << " "
                 << nodeIndexFromIterator(eIt, 0) << " " << nodeIndexFromIterator(eIt, 1) << " "
                 << nodeIndexFromIterator(eIt, 3) << " " << nodeIndexFromIterator(eIt, 2) << " "
                 << nodeIndexFromIterator(eIt, 4) << " " << nodeIndexFromIterator(eIt, 5) << " "
                 << nodeIndexFromIterator(eIt, 7) << " " << nodeIndexFromIterator(eIt, 6);
          else if (7 == element_type)
            file << " "
                 << nodeIndexFromIterator(eIt, 0) << " " << nodeIndexFromIterator(eIt, 1) << " "
                 << nodeIndexFromIterator(eIt, 3) << " " << nodeIndexFromIterator(eIt, 2) << " "
                 << nodeIndexFromIterator(eIt, 4);
          else {
            for (int k = 0; k < eIt->geometry().corners(); ++k)
              file << " " << nodeIndexFromIterator(eIt, k);
          }

          file << std::endl;

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
      VertexIterator vIt    = gv.template begin<dim>();
      VertexIterator vEndIt = gv.template end<dim>();

      for (; vIt != vEndIt; ++vIt) {
        typename VertexIterator::Entity::Geometry::GlobalCoordinate globalCoord = vIt->geometry().center();
        int nodeIndex = gv.indexSet().index(*vIt)+1; // Start counting indices by "1".

        if (1 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << 0 << " " << 0 << std::endl;
        else if (2 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << globalCoord[1] << " " << 0 << std::endl;
        else // (3 == dimWorld)
          file << nodeIndex << " " << globalCoord[0] << " " << globalCoord[1] << " " << globalCoord[2] << std::endl;
      }
    }


  public:
    /** \brief Constructor expecting GridView of Grid to be written.
        \param gridView GridView that will be used in write(const std::string&).
    */
    GmshWriter(const GridView& gridView) : gv(gridView) {}

    /** \brief Write given grid in Gmsh 2.0 compatible ASCII file.
        \param fileName Path of file. write(const std::string&) does not attach a ".msh"-extension by itself.

        Opens the file with given name and path, stores the element data of the grid
        and closes the file when done.

        Boundary-Segments of the grid are ignored.

        Throws an IOError if file could not be opened or an unsupported element type is
        encountered.
    */
    void write(const std::string& fileName) const {
      std::ofstream file(fileName.c_str());

      if (!file.is_open())
        DUNE_THROW(Dune::IOError, "Could not open " << fileName << " with write access.");

      // Output Header
      file << "$MeshFormat" << std::endl
           << "2.0 0 " << sizeof(double) << std::endl // "2.0" for "version 2.0", "0" for ASCII
           << "$EndMeshFormat" << std::endl;


      // Output Nodes
      const size_t number_of_nodes = gv.size(dim);
      file << "$Nodes" << std::endl
           << number_of_nodes << std::endl;

      outputNodes(file);

      file << "$EndNodes" << std::endl;


      // Output Elements
      const size_t number_of_elements = gv.size(0);
      file << "$Elements" << std::endl
           << number_of_elements << std::endl;

      outputElements(file);

      file << "$EndElements" << std::endl; // Add an additional new line for good measure.


      file.close();
    }
  };

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_GMSHWRITER_HH
