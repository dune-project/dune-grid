// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <algorithm>
#include <memory>
#include <vector>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/common/stdstreams.hh>

#include <amiramesh/AmiraMesh.h>

#if HAVE_PSURFACE
#include "psurface/AmiraMeshIO.h"
#endif



template <class GridType>
template <class DiscFuncType>
void Dune::AmiraMeshReader<GridType>::readFunction(DiscFuncType& f, const std::string& filename)
{
  // f may be a block vector
  const int blocksize = DiscFuncType::block_type::dimension;

  // /////////////////////////////////////////////////////
  // Load the AmiraMesh file
  std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));

  if(!am)
    DUNE_THROW(IOError, "Could not open AmiraMesh file: " << filename);

  int i, j;

  f = 0;

  bool datafound=false;

  // We allow P1 fields defined on the whole grid and fields defined
  // only on the boundary.  We now check the file and proceed accordingly
  if (!am->findData("Nodes", HxFLOAT, blocksize, "Data") &&
      !am->findData("Nodes", HxDOUBLE, blocksize, "Data")) {

    // get the data field
    AmiraMesh::Data* am_ValueData =  am->findData("Nodes", HxFLOAT, blocksize, "values");
    if (am_ValueData) {
      datafound = true;

      if ((int) f.size() < am->nElements("Nodes"))
        DUNE_THROW(IOError, "When reading data from a surface field the "
                   << "array you provide has to have at least the size of the surface!");

      float* am_values_float = (float*) am_ValueData->dataPtr();

      for (i=0; i<am->nElements("Nodes"); i++) {
        for (j=0; j<blocksize; j++)
          f[i][j] = am_values_float[i*blocksize+j];

      }

    } else {
      am_ValueData =  am->findData("Nodes", HxDOUBLE, blocksize, "values");
      if (am_ValueData) {
        datafound = true;

        if ((int) f.size() < am->nElements("Nodes"))
          DUNE_THROW(IOError, "When reading data from a surface field your the "
                     << "array you provide has to have at least the size of the surface!");

        for (i=0; i<blocksize*am->nElements("Nodes"); i++) {
          for (j=0; j<blocksize; j++)
            f[i][j] = ((double*)am_ValueData->dataPtr())[i*blocksize+j];
        }

      }
    }

  } else {

    // get the data field
    AmiraMesh::Data* am_ValueData =  am->findData("Nodes", HxFLOAT, blocksize, "Data");
    if (am_ValueData) {
      datafound = true;

      float* am_values_float = (float*) am_ValueData->dataPtr();
      f.resize(am->nElements("Nodes"));

      for (i=0; i<am->nElements("Nodes"); i++)
        for (j=0; j<blocksize; j++)
          f[i][j] = am_values_float[i*blocksize+j];

    } else {
      am_ValueData =  am->findData("Nodes", HxDOUBLE, blocksize, "Data");
      if (am_ValueData) {
        datafound = true;
        f.resize(am->nElements("Nodes"));

        for (i=0; i<am->nElements("Nodes"); i++)
          for (j=0; j<blocksize; j++)
            f[i][j] = ((double*)am_ValueData->dataPtr())[i*blocksize+j];
      } else
        DUNE_THROW(IOError, "No data found in the file!");
    }
  }
  // No proper P1 function has been found for the grid maybe it's P0 ?
  if (!datafound and (am->findData("Triangles", HxFLOAT, blocksize, "Data") or am->findData("Triangles", HxDOUBLE, blocksize, "Data")))
  {
    // get the data field
    AmiraMesh::Data* am_ValueData =  am->findData("Triangles", HxFLOAT, blocksize, "Data");
    if (am_ValueData)
    {
      datafound = true;

      float* am_values_float = (float*) am_ValueData->dataPtr();
      f.resize(am->nElements("Triangles"));

      for (i=0; i<am->nElements("Triangles"); i++)
        for (j=0; j<blocksize; j++)
          f[i][j] = am_values_float[i*blocksize+j];

    }
    else
    {
      am_ValueData =  am->findData("Triangles", HxDOUBLE, blocksize, "Data");
      if (am_ValueData)
      {
        datafound = true;
        f.resize(am->nElements("Triangles"));

        for (i=0; i<am->nElements("Triangles"); i++)
          for (j=0; j<blocksize; j++)
            f[i][j] = ((double*)am_ValueData->dataPtr())[i*blocksize+j];
      }
      else
        DUNE_THROW(IOError, "No data found in the file!");
    }
  }
  if (!datafound)
    DUNE_THROW(IOError, "No data found in the file!");

  dverb << "Data field " << filename << " loaded successfully!" << std::endl;

}


// Create the domain from an explicitly given boundary description
template <class GridType>
void Dune::AmiraMeshReader<GridType>::createDomain(GridFactory<GridType>& factory,
                                                   const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary)
{
#if HAVE_PSURFACE
  if (dim!=3)
    DUNE_THROW(NotImplemented, "AmiraMeshReader with PSurface support only for 3d grids");

  psurface::PSurface<2, ctype>* psurface = reinterpret_cast<psurface::PSurface<2, ctype>* >(boundary->getPSurfaceObject());

  if (!psurface->hasUpToDatePointLocationStructure)
    psurface->createPointLocationStructure();

  int noOfSegments = psurface->getNumTriangles();
  if(noOfSegments <= 0)
    DUNE_THROW(IOError, "no segments found");

  int noOfNodes = psurface->getNumVertices();
  if(noOfNodes <= 0)
    DUNE_THROW(IOError, "No nodes found");

  const int dim = GridType::dimension;

  for(int i = 0; i < noOfSegments; i++) {

    // Gets the vertices of a boundary segment
    std::vector<unsigned int> vertices(3);
    vertices[0] = psurface->triangles(i).vertices[0];
    vertices[1] = psurface->triangles(i).vertices[1];
    vertices[2] = psurface->triangles(i).vertices[2];

    factory.insertBoundarySegment(vertices,
                                  std::shared_ptr<BoundarySegment<dim,dim> >(new typename PSurfaceBoundary<dim-1, ctype>::PSurfaceBoundarySegment(boundary,i)));

  }

  Dune::dinfo << noOfSegments << " segments from PSurfaceBoundary object created!" << std::endl;

#endif // #define HAVE_PSURFACE
}


template <class GridType>
Dune::ToUniquePtr<GridType> Dune::AmiraMeshReader<GridType>::read(const std::string& filename,
                                                const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary)
{
#if ! HAVE_PSURFACE
  DUNE_THROW(IOError, "Dune has not been built with support for the "
             << " psurface library!");
#else
  dverb << "This is the AmiraMesh reader for " << className<GridType>() << std::endl;

  // Create a grid factory
  GridFactory<GridType> factory;

  // /////////////////////////////////////////////////////
  // Load the AmiraMesh file
  // /////////////////////////////////////////////////////
  std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));

  if(!am)
    DUNE_THROW(IOError, "Could not open AmiraMesh file " << filename);

  if (am->findData("Hexahedra", HxINT32, 8, "Nodes")) {

    // Load a domain from an AmiraMesh hexagrid file
    Dune::dwarn << "Hexahedral grids with a parametrized boundary are not supported!" << std::endl;
    Dune::dwarn << "I will therefore ignore the boundary parametrization." << std::endl;

  } else {

    // Load domain from an AmiraMesh tetragrid file
    createDomain(factory, boundary);

  }

  // read and build the grid
  buildGrid(factory, am.get());

  return factory.createGrid();
#endif // #define HAVE_PSURFACE
}


template <class GridType>
void Dune::AmiraMeshReader<GridType>::read(GridType& grid,
                                           const std::string& filename,
                                           const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary)
{
#if ! HAVE_PSURFACE
  DUNE_THROW(IOError, "Dune has not been built with support for the "
             << " psurface library!");
#else
  dverb << "This is the AmiraMesh reader for" << className<GridType>() << "!" << std::endl;

  // Create a grid factory
  GridFactory<GridType> factory(&grid);

  // /////////////////////////////////////////////////////
  // Load the AmiraMesh file
  // /////////////////////////////////////////////////////
  std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));

  if(!am)
    DUNE_THROW(IOError, "Could not open AmiraMesh file " << filename);

  if (am->findData("Hexahedra", HxINT32, 8, "Nodes")) {

    // Load a domain from an AmiraMesh hexagrid file
    Dune::dwarn << "Hexahedral grids with a parametrized boundary are not supported!" << std::endl;
    Dune::dwarn << "I will therefore ignore the boundary parametrization." << std::endl;

  } else {

    // Load domain from an AmiraMesh tetragrid file
    createDomain(factory, boundary);

  }

  // read and build the grid
  buildGrid(factory, am.get());

  factory.createGrid();
#endif // #define HAVE_PSURFACE
}


template <class GridType>
Dune::ToUniquePtr<GridType> Dune::AmiraMeshReader<GridType>::read(const std::string& filename)
{
  static const int dim      = GridType::dimension;
  static const int dimworld = GridType::dimensionworld;

  static_assert(dim==dimworld, "AmiraMesh can only be read for grids with dim==dimworld!");

  // Create a grid factory
  GridFactory<GridType> factory;

  // Load the AmiraMesh file
  std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));

  if(!am)
    DUNE_THROW(IOError, "read: Could not open AmiraMesh file " << filename);

  // ////////////////////////////////////////////////////
  //   Build the grids
  // ////////////////////////////////////////////////////
  if (dim==3) {

    buildGrid(factory, am.get());

  } else {

    build2dGrid(factory, am.get());

  }

  return factory.createGrid();
}

template <class GridType>
void Dune::AmiraMeshReader<GridType>::read(GridType& grid,
                                           const std::string& filename)
{
  static const int dim      = GridType::dimension;
  static const int dimworld = GridType::dimensionworld;

  if (dim!=dimworld)
    DUNE_THROW(IOError, "AmiraMesh can only be read for grids with dim==dimworld!");

  dverb << "This is the AmiraMesh reader for file '" << filename << "'!" << std::endl;

  // Create a grid factory
  GridFactory<GridType> factory(&grid);

  // Load the AmiraMesh file
  std::unique_ptr<AmiraMesh> am(AmiraMesh::read(filename.c_str()));
  if(!am)
    DUNE_THROW(IOError, "read: Could not open AmiraMesh file " << filename);

  // ////////////////////////////////////////////////////
  //   Build the grids
  // ////////////////////////////////////////////////////
  if (dim==3) {

    buildGrid(factory, am.get());

  } else {

    build2dGrid(factory, am.get());

  }

  factory.createGrid();
}

template <class GridType>
void Dune::AmiraMeshReader<GridType>::build2dGrid(GridFactory<GridType>& factory, AmiraMesh* am)
{
  static const int dimworld = GridType::dimensionworld;

  // ////////////////////////////////////////////////////
  //   Here we now the grid must be 2d
  // ////////////////////////////////////////////////////

  // Determine whether grid contains only triangles
  bool containsOnlyTriangles = am->findData("Triangles", HxINT32, 3, "Nodes");

  // get the different data fields
  AmiraMesh::Data* am_coordinateData =  am->findData("Nodes", HxFLOAT, 2, "Coordinates");
  if (!am_coordinateData)
    DUNE_THROW(IOError, "2D AmiraMesh loader: field 'Nodes' not found!");

  float* am_node_coordinates = (float*) am_coordinateData->dataPtr();

  // Get the element list
  int*  elemData = 0;

  if (containsOnlyTriangles) {
    AmiraMesh::Data* triangleData = am->findData("Triangles", HxINT32, 3, "Nodes");
    if (triangleData)
      elemData = (int*)triangleData->dataPtr();
    else
      DUNE_THROW(IOError, "2D AmiraMesh loader: field 'Triangles' not found!");

  } else {
    AmiraMesh::Data* elementData = am->findData("Quadrilaterals", HxINT32, 4, "Nodes");
    if (elementData) {
      elemData = (int*)elementData->dataPtr();
    } else
      DUNE_THROW(IOError, "2D AmiraMesh loader: field 'Quadrilaterals' not found!");
  }


  int noOfNodes = am->nElements("Nodes");
  int noOfElem  = (containsOnlyTriangles) ? am->nElements("Triangles") : am->nElements("Quadrilaterals");

  Dune::dinfo << "AmiraMesh contains " << noOfNodes << " nodes and "
              << noOfElem << " elements\n";

  // Insert interior nodes
  for(int i=0; i < noOfNodes; i++) {

    FieldVector<double,dimworld> nodePos;
    nodePos[0] = am_node_coordinates[2*i];
    nodePos[1] = am_node_coordinates[2*i+1];

    factory.insertVertex(nodePos);

  }

  int noOfCreatedElem = 0;
  for (int i=0; i < noOfElem; i++) {

    if (containsOnlyTriangles) {

      std::vector<unsigned int> cornerIDs(3);

      /* only triangles */
      cornerIDs[0] = elemData[3*i]-1;
      cornerIDs[1] = elemData[3*i+1]-1;
      cornerIDs[2] = elemData[3*i+2]-1;

      factory.insertElement(GeometryTypes::triangle, cornerIDs);

    } else {

      if (elemData[4*i+2]==elemData[4*i+3]) {
        // Triangle within a quadrilateral grid file
        std::vector<unsigned int> cornerIDs(3);

        /* a triangle in a quadrilateral file */
        cornerIDs[0] = elemData[4*i]-1;
        cornerIDs[1] = elemData[4*i+1]-1;
        cornerIDs[2] = elemData[4*i+2]-1;

        factory.insertElement(GeometryTypes::triangle, cornerIDs);

      } else {

        std::vector<unsigned int> cornerIDs(4);

        /* a true quadrilateral */
        cornerIDs[0] = elemData[4*i]-1;
        cornerIDs[1] = elemData[4*i+1]-1;
        cornerIDs[2] = elemData[4*i+3]-1;
        cornerIDs[3] = elemData[4*i+2]-1;

        factory.insertElement(GeometryTypes::quadrilateral, cornerIDs);

      }

    }

    noOfCreatedElem++;

  }

  Dune::dinfo << "amiraloadmesh: " << noOfCreatedElem << " elements created" << std::endl;

}

template <class GridType>
void Dune::AmiraMeshReader<GridType>::buildGrid(Dune::GridFactory<GridType>& factory,
                                                AmiraMesh* am)
{
  static const int dimworld = GridType::dimensionworld;

  bool isTetraGrid = am->findData("Tetrahedra", HxINT32, 4, "Nodes");

  float* am_node_coordinates_float = NULL;
  double* am_node_coordinates_double = NULL;

  // get the different data fields
  AmiraMesh::Data* am_coordinateData =  am->findData("Nodes", HxFLOAT, 3, "Coordinates");
  if (am_coordinateData)
    am_node_coordinates_float = (float*) am_coordinateData->dataPtr();
  else {
    am_coordinateData =  am->findData("Nodes", HxDOUBLE, 3, "Coordinates");
    if (am_coordinateData)
      am_node_coordinates_double = (double*) am_coordinateData->dataPtr();
    else
      DUNE_THROW(IOError, "No vertex coordinates found in the file!");

  }


  AmiraMesh::Data* elementData = (isTetraGrid)
                                 ? am->findData("Tetrahedra", HxINT32, 4, "Nodes")
                                 : am->findData("Hexahedra", HxINT32, 8, "Nodes");

  int*  elemData = (int*)elementData->dataPtr();

  int noOfNodes = am->nElements("Nodes");

  Dune::dinfo << "AmiraMesh has " << noOfNodes << " total nodes." << std::endl;

  int noOfElem = (isTetraGrid)
                 ? am->nElements("Tetrahedra")
                 : am->nElements("Hexahedra");

  // //////////////////////////////////////
  //   Insert interior nodes
  // //////////////////////////////////////
  assert(am_node_coordinates_float || am_node_coordinates_double);
  for(int i=0; i < noOfNodes; i++) {

    FieldVector<double,dimworld> nodePos;

    if (am_node_coordinates_float) {
      nodePos[0] = am_node_coordinates_float[3*i];
      nodePos[1] = am_node_coordinates_float[3*i+1];
      nodePos[2] = am_node_coordinates_float[3*i+2];
    } else {
      nodePos[0] = am_node_coordinates_double[3*i];
      nodePos[1] = am_node_coordinates_double[3*i+1];
      nodePos[2] = am_node_coordinates_double[3*i+2];
    }

    factory.insertVertex(nodePos);

  }

  /* all inner nodes are inserted , now we insert the elements */
  for(int i=0; i < noOfElem; i++) {

    const int* thisElem = elemData + (i* ((isTetraGrid) ? 4 : 8));

    if (isTetraGrid) {

      int numberOfCorners = 4;
      std::vector<unsigned int> cornerIDs(numberOfCorners);

      for (int j=0; j<numberOfCorners; j++)
        cornerIDs[j] = elemData[numberOfCorners*i+j]-1;

      factory.insertElement(GeometryTypes::tetrahedron, cornerIDs);

    } else {

      if (thisElem[2]==thisElem[3]
          && thisElem[4]==thisElem[5]
          && thisElem[5]==thisElem[6]
          && thisElem[6]==thisElem[7]) {

        // Tetrahedron
        std::vector<unsigned int> cornerIDs(4);

        cornerIDs[0] = thisElem[0]-1;
        cornerIDs[1] = thisElem[1]-1;
        cornerIDs[2] = thisElem[2]-1;
        cornerIDs[3] = thisElem[4]-1;

        factory.insertElement(GeometryTypes::tetrahedron, cornerIDs);

      }else if (thisElem[4]==thisElem[5] && thisElem[5]==thisElem[6]
                && thisElem[6]==thisElem[7]) {

        // Pyramid
        std::vector<unsigned int> cornerIDs(5);

        cornerIDs[0] = thisElem[0]-1;
        cornerIDs[1] = thisElem[1]-1;
        cornerIDs[2] = thisElem[2]-1;
        cornerIDs[3] = thisElem[3]-1;
        cornerIDs[4] = thisElem[4]-1;

        factory.insertElement(GeometryTypes::pyramid, cornerIDs);

      } else if (thisElem[1]==thisElem[2] && thisElem[5]==thisElem[6]) {

        // Prism
        std::vector<unsigned int> cornerIDs(6);

        cornerIDs[0] = thisElem[0]-1;
        cornerIDs[1] = thisElem[1]-1;
        cornerIDs[2] = thisElem[3]-1;
        cornerIDs[3] = thisElem[4]-1;
        cornerIDs[4] = thisElem[5]-1;
        cornerIDs[5] = thisElem[7]-1;

        factory.insertElement(GeometryTypes::prism, cornerIDs);

      } else if (thisElem[2]==thisElem[3] && thisElem[6]==thisElem[7]) {

        std::vector<unsigned int> cornerIDs(6);

        cornerIDs[0] = thisElem[0]-1;
        cornerIDs[1] = thisElem[1]-1;
        cornerIDs[2] = thisElem[2]-1;
        cornerIDs[3] = thisElem[4]-1;
        cornerIDs[4] = thisElem[5]-1;
        cornerIDs[5] = thisElem[6]-1;

        factory.insertElement(GeometryTypes::prism, cornerIDs);

      } else {

        int numberOfCorners = 8;
        std::vector<unsigned int> cornerIDs(numberOfCorners);

        cornerIDs[0] = elemData[numberOfCorners*i+0]-1;
        cornerIDs[1] = elemData[numberOfCorners*i+1]-1;
        cornerIDs[2] = elemData[numberOfCorners*i+3]-1;
        cornerIDs[3] = elemData[numberOfCorners*i+2]-1;
        cornerIDs[4] = elemData[numberOfCorners*i+4]-1;
        cornerIDs[5] = elemData[numberOfCorners*i+5]-1;
        cornerIDs[6] = elemData[numberOfCorners*i+7]-1;
        cornerIDs[7] = elemData[numberOfCorners*i+6]-1;

        factory.insertElement(GeometryTypes::hexahedron, cornerIDs);

      }

    }

  }

  Dune::dinfo << "AmiraMesh reader: " << noOfElem << " elements created.\n";

}
