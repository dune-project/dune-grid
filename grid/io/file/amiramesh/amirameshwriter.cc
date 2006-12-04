// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <amiramesh/AmiraMesh.h>

#include <algorithm>

// ///////////////////////////////////////////////////////////
//   Write the leaf level of the grid as an AmiraMesh file
// ///////////////////////////////////////////////////////////
template<class GridType, class IndexSetType>
void Dune::AmiraMeshWriter<GridType,IndexSetType>::writeGrid(const GridType& grid,
                                                             const std::string& filename,
                                                             const IndexSetType& indexSet)
{
  if ((dim!=2 && dim!=3) || int(dim) != int(GridType::dimensionworld))
    DUNE_THROW(IOError, "You can only write grids as AmiraMesh if dim==dimworld==2"
               << " or dim==dimworld==3.");

  // Find out whether the grid contains only tetrahedra.  If yes, then
  // it is written in TetraGrid format.  If not, it is written in
  // hexagrid format.
  bool containsOnlySimplices =
    (indexSet.geomTypes(0).size()==1)
    && (indexSet.geomTypes(0)[0].isSimplex());

  int maxVerticesPerElement = (dim==3)
                              ? ((containsOnlySimplices) ? 4 : 8)
                              : ((containsOnlySimplices) ? 3 : 4);

  int noOfNodes = indexSet.size(dim);
  int noOfElem  = indexSet.size(0);

  // create amiramesh object
  AmiraMesh am;

  // Set the appropriate content type
  if (dim==2)
    am.parameters.set("ContentType", "HxTriangularGrid");

  // write grid vertex coordinates
  AmiraMesh::Location* geo_nodes = new AmiraMesh::Location("Nodes", noOfNodes);
  am.insert(geo_nodes);

  AmiraMesh::Data* geo_node_data = new AmiraMesh::Data("Coordinates", geo_nodes,
                                                       McPrimType::mc_float, dim);
  am.insert(geo_node_data);

  typedef typename IndexSetType::template Codim<dim>::template Partition<All_Partition>::Iterator VertexIterator;
  VertexIterator vertex    = indexSet.template begin<dim,All_Partition>();
  VertexIterator endvertex = indexSet.template end<dim,All_Partition>();

  for (; vertex!=endvertex; ++vertex) {

    int index = indexSet.index (*vertex);
    const FieldVector<double, dim>& coords = vertex->geometry()[0];

    // Copy coordinates
    for (int i=0; i<dim; i++)
      ((float*)geo_node_data->dataPtr())[dim*index+i] = coords[i];

  }

  /* write element section to file */
  AmiraMesh::Location* element_loc = NULL;

  if (dim==3) {

    if (containsOnlySimplices)
      element_loc = new AmiraMesh::Location("Tetrahedra", noOfElem);
    else
      element_loc = new AmiraMesh::Location("Hexahedra", noOfElem);

  } else {

    if (containsOnlySimplices)
      element_loc = new AmiraMesh::Location("Triangles", noOfElem);
    else
      element_loc = new AmiraMesh::Location("Quadrilaterals", noOfElem);

  }

  am.insert(element_loc);

  AmiraMesh::Data* element_data = new AmiraMesh::Data("Nodes", element_loc,
                                                      McPrimType::mc_int32, maxVerticesPerElement);
  am.insert(element_data);

  int *dPtr = (int*)element_data->dataPtr();

  typedef typename IndexSetType::template Codim<0>::template Partition<All_Partition>::Iterator ElementIterator;
  ElementIterator eIt    = indexSet.template begin<0, All_Partition>();
  ElementIterator eEndIt = indexSet.template end<0, All_Partition>();

  if (dim==3) {

    // //////////////////////////////////////////////////
    //   Write elements of a 3D-grid
    // //////////////////////////////////////////////////

    if (containsOnlySimplices) {

      for (int i=0; eIt!=eEndIt; ++eIt, i++) {

        for (int j=0; j<4; j++)
          dPtr[i*4+j] = indexSet.template subIndex<dim>(*eIt,j)+1;

      }

    } else {

      for (int i=0; eIt!=eEndIt; ++eIt, i++) {

        GeometryType type = eIt->geometry().type();

        if (type.isHexahedron()) {

          const int hexaReordering[8] = {0, 1, 3, 2, 4, 5, 7, 6};
          for (int j=0; j<8; j++)
            dPtr[8*i + j] = indexSet.template subIndex<dim>(*eIt, hexaReordering[j])+1;

        } else if (type.isPrism()) {

          const int prismReordering[8] = {0, 1, 1, 2, 3, 4, 4, 5};
          for (int j=0; j<8; j++)
            dPtr[8*i + j] = indexSet.template subIndex<dim>(*eIt, prismReordering[j])+1;

        } else if (type.isPyramid()) {

          const int pyramidReordering[8] = {0, 1, 2, 3, 4, 4, 4, 4};
          for (int j=0; j<8; j++)
            dPtr[8*i + j] = indexSet.template subIndex<dim>(*eIt, pyramidReordering[j])+1;

        } else if (type.isTetrahedron()) {

          const int tetraReordering[8] = {0, 1, 2, 2, 3, 3, 3, 3};
          for (int j=0; j<8; j++)
            dPtr[8*i + j] = indexSet.template subIndex<dim>(*eIt, tetraReordering[j])+1;

        } else {

          DUNE_THROW(NotImplemented, "Unknown element type encountered");

        }

      }

    }

  } else {

    for (int i=0; eIt!=eEndIt; ++eIt, i++) {

      GeometryType type = eIt->geometry().type();

      if (type.isQuadrilateral()) {

        dPtr[i*4+0] = indexSet.template subIndex<dim>(*eIt, 0)+1;
        dPtr[i*4+1] = indexSet.template subIndex<dim>(*eIt, 1)+1;
        dPtr[i*4+2] = indexSet.template subIndex<dim>(*eIt, 3)+1;
        dPtr[i*4+3] = indexSet.template subIndex<dim>(*eIt, 2)+1;

      } else if (type.isTriangle()) {

        for (int j=0; j<3; j++)
          dPtr[i*maxVerticesPerElement+j] = indexSet.template subIndex<dim>(*eIt, j)+1;

        // If 4 vertices are expected per element use the last value
        // to fill up the remaining slots
        if (maxVerticesPerElement==4)
          dPtr[i*4+3] = dPtr[i*4+2];

      } else {

        DUNE_THROW(IOError, "Elements of type " << type
                                                << " cannot be written to 2d AmiraMesh files!");
      }

    }

  }

  // write material section to grid file
  AmiraMesh::Data* element_materials = new AmiraMesh::Data("Materials", element_loc, McPrimType::mc_uint8, 1);
  am.insert(element_materials);

  for(int i=0; i<noOfElem; i++)
    ((unsigned char*)element_materials->dataPtr())[i] = 0;

  // Actually write the file
  if(!am.write(filename.c_str(), 1))
    DUNE_THROW(IOError, "Writing geometry file failed!");

  std::cout << "Grid written successfully to: " << filename << std::endl;
}

template<class GridType, class IndexSetType>
template<class VectorType>
void Dune::AmiraMeshWriter<GridType,IndexSetType>::writeBlockVector(const VectorType& f,
                                                                    const std::string& filename,
                                                                    const IndexSetType& indexSet)
{
  // Find out whether the grid contains only tetrahedra.  If yes, then
  // it is written in TetraGrid format.  If not, it is written in
  // hexagrid format.
  bool containsOnlyTetrahedra =
    (indexSet.geomTypes(0).size()==1)
    && (indexSet.geomTypes(0)[0].isSimplex());


  // Get number of components
  const int ncomp = VectorType::block_type::size;

  // Create AmiraMesh object
  AmiraMesh am;

  // Set the appropriate content type for 2D grid data
  if (GridType::dimension==2)
    am.parameters.set("ContentType", "HxTriangularData");

  if (!containsOnlyTetrahedra && dim==3) {

    AmiraMesh::Location* hexa_loc = new AmiraMesh::Location("Hexahedra", indexSet.size(0));
    am.insert(hexa_loc);

  }

  AmiraMesh::Location* amLocation;
  if (f.size()==indexSet.size(dim)) {

    // P1 data
    amLocation = new AmiraMesh::Location("Nodes", f.size());

  } else if (f.size()==indexSet.size(0)) {

    // P0 data
    if (containsOnlyTetrahedra)
      amLocation = new AmiraMesh::Location((dim==2) ? "Triangles" : "Tetrahedra", f.size());
    else
      amLocation = new AmiraMesh::Location("Hexahedra", f.size());

  } else
    DUNE_THROW(IOError, "BlockVector doesn't match the grid!");

  am.insert(amLocation);

  AmiraMesh::Data* nodeData = new AmiraMesh::Data("Data", amLocation, McPrimType::mc_double, ncomp);
  am.insert(nodeData);


  AmiraMesh::Field* nodeField;

  if (f.size()==indexSet.size(0)) {
    // P0 data
    nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double,
                                     AmiraMesh::t_constant, nodeData);

  } else if (containsOnlyTetrahedra || dim==2) {
    nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double,
                                     AmiraMesh::t_linear, nodeData);
  } else {

    nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double,
                                     AmiraMesh::t_trilinear, nodeData);
  }

  am.insert(nodeField);


  // write the data into the AmiraMesh object
  typedef typename VectorType::ConstIterator Iterator;
  Iterator dit    = f.begin();
  Iterator ditend = f.end();

  int i=0;
  for (; dit!=ditend; ++dit) {

    for (int j=0; j<ncomp; j++)
      ((double*)nodeData->dataPtr())[i++] = (*dit)[j];

  }

  // actually save the solution file
  // (the 1 means: ascii)
  if (!am.write(filename.c_str(), 1) )
    DUNE_THROW(IOError, "An error has occured writing file " << filename);

  std::cout << "Solution written successfully to: " << filename << std::endl;
}
