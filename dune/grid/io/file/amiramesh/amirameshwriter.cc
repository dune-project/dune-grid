// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <amiramesh/AmiraMesh.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <dune/grid/common/virtualrefinement.hh>

template<class GridView>
void Dune::AmiraMeshWriter<GridView>::addGrid(const GridView& gridView,
                                              bool splitAll)
{

  typedef typename GridView::Grid::ctype ct;
  typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
  typedef typename GridView::template Codim<0>::Iterator ElementIterator;
  typedef Dune::VirtualRefinement<dim, ct> Refinement;
  typedef typename Refinement::VertexIterator vIterator;
  typedef typename Refinement::ElementIterator eIterator;
  typedef typename Refinement::VertexIterator::CoordVector Coordinate;
  typedef typename Refinement::ElementIterator::IndexVector IndexVector;

  const typename GridView::IndexSet& indexSet = gridView.indexSet();

  if ((dim!=2 && dim!=3) || int(dim) != int(GridView::dimensionworld))
    DUNE_THROW(IOError, "You can only write grids as AmiraMesh if dim==dimworld==2"
               << " or dim==dimworld==3.");

  // Set the appropriate content type
  if (dim==2)
    amiramesh_.parameters.set("ContentType", "HxTriangularGrid");

  // ///////////////////////////////////////////
  //   Write grid vertices
  // ///////////////////////////////////////////
  int noOfNodes = indexSet.size(dim);
  int noOfElements;

  AmiraMesh::Location* geo_nodes = new AmiraMesh::Location("Nodes", noOfNodes);
  amiramesh_.insert(geo_nodes);

  AmiraMesh::Data* geo_node_data = new AmiraMesh::Data("Coordinates", geo_nodes,
                                                       McPrimType::mc_float, dim);
  amiramesh_.insert(geo_node_data);

  VertexIterator vertex    = gridView.template begin<dim>();
  VertexIterator endvertex = gridView.template end<dim>();

  //needed later to compare indices of refinement and indexset
  std::vector<Coordinate> vertices_coords(noOfNodes);


  for (; vertex!=endvertex; ++vertex) {

    int index = indexSet.template index<dim>(*vertex);

    vertices_coords[index]=vertex->geometry().corner(0);

    // Copy coordinates
    for (int i=0; i<dim; i++)
      ((float*)geo_node_data->dataPtr())[dim*index+i] = vertex->geometry().corner(0)[i];

  }


  /* write element section to file */
  AmiraMesh::Location* elementLocation = NULL;

  // ////////////////////////////////////////////////////////////////////
  //   Split up all elements into simplices, if requested, because
  //   Amira doesn't support all kind of grids.
  // ////////////////////////////////////////////////////////////////////



  if (splitAll) {
    Dune::GeometryType coerceTo(Dune::GeometryType::simplex,dim);
    noOfElements = 0;
    int count=0;

    for (size_t i=0; i<indexSet.geomTypes(0).size(); i++) {

      if (indexSet.geomTypes(0)[i].isSimplex())
        count=1;
      else {
        Refinement & refinement = Dune::buildRefinement<dim, ct>(indexSet.geomTypes(0)[i],coerceTo);
        count = refinement.nElements(0);
      }

      noOfElements += count * indexSet.size(indexSet.geomTypes(0)[i]);

    }

    // write element section to file

    int VerticesPerElement = dim + 1;

    if (dim==2)
      elementLocation = new AmiraMesh::Location("Triangles", noOfElements);
    else
      elementLocation = new AmiraMesh::Location("Tetrahedra", noOfElements);

    amiramesh_.insert(elementLocation);

    AmiraMesh::Data* element_data = new AmiraMesh::Data("Nodes", elementLocation, McPrimType::mc_int32, VerticesPerElement);
    amiramesh_.insert(element_data);

    int *dPtr = (int*)element_data->dataPtr();

    ElementIterator eIt    = gridView.template begin<0>();
    ElementIterator eEndIt = gridView.template end<0>();

    for (int i=0; eIt!=eEndIt; ++eIt) {
      if (eIt->type().isSimplex()) {
        for (int j=0; j<VerticesPerElement; j++)
          // The +1 is added because AmiraMesh numbers vertices starting from 1
          dPtr[i*VerticesPerElement+j] = indexSet.subIndex(*eIt,j,dim)+1;
        i++;
      }
      else {

        Refinement & refinement = Dune::buildRefinement<dim, ct>(eIt->type(),coerceTo);

        eIterator eSubEnd = refinement.eEnd(0);
        eIterator eSubIt = refinement.eBegin(0);
        IndexVector vertexIds;

        //Have to do this, because Refinement indices of corners don't match indexSet indices of the corners of that entity
        //So we have to check equality of two indices by checking equality of the coordinates of the corners

        //coordinates of nodes using refinement indices
        std::vector<Coordinate> vertices(refinement.nVertices(0));

        vIterator vEnd = refinement.vEnd(0);
        for(vIterator vIt = refinement.vBegin(0); vIt != vEnd; ++vIt)
          vertices[vIt.index()]=(eIt->geometry().global(vIt.coords()));

        for( ; eSubIt != eSubEnd; ++eSubIt) {

          vertexIds = eSubIt.vertexIndices();

          for (int j=0; j<VerticesPerElement; j++) {

            //better way to do the comparison ?
            for (int m=0; m<noOfNodes; m++) {
              if (vertices[vertexIds[j]]==vertices_coords[m]) {
                dPtr[i*VerticesPerElement+j]=m+1;
                break;
              }
            }

          }
          i++;

        }

      }

    }

  } else {
    // Find out whether the grid contains only tetrahedra.  If yes, then
    // it is written in TetraGrid format.  If not, it is written in
    // hexagrid format.
    bool containsOnlySimplices =
      (indexSet.geomTypes(0).size()==1)
      && (indexSet.geomTypes(0)[0].isSimplex());

    int maxVerticesPerElement = (dim==3)
                                ? ((containsOnlySimplices) ? 4 : 8)
                                : ((containsOnlySimplices) ? 3 : 4);

    noOfElements  = indexSet.size(0);

    // write element section to file
    if (dim==3) {

      if (containsOnlySimplices)
        elementLocation = new AmiraMesh::Location("Tetrahedra", noOfElements);
      else
        elementLocation = new AmiraMesh::Location("Hexahedra", noOfElements);

    } else {

      if (containsOnlySimplices)
        elementLocation = new AmiraMesh::Location("Triangles", noOfElements);
      else
        elementLocation = new AmiraMesh::Location("Quadrilaterals", noOfElements);

    }

    amiramesh_.insert(elementLocation);

    AmiraMesh::Data* element_data = new AmiraMesh::Data("Nodes", elementLocation,
                                                        McPrimType::mc_int32, maxVerticesPerElement);
    amiramesh_.insert(element_data);

    int *dPtr = (int*)element_data->dataPtr();

    ElementIterator eIt    = gridView.template begin<0>();
    ElementIterator eEndIt = gridView.template end<0>();

    if (dim==3) {

      // //////////////////////////////////////////////////
      //   Write elements of a 3D-grid
      // //////////////////////////////////////////////////

      if (containsOnlySimplices) {

        for (int i=0; eIt!=eEndIt; ++eIt, i++) {

          for (int j=0; j<4; j++)
            dPtr[i*4+j] = indexSet.subIndex(*eIt,j,dim)+1;

        }

      } else {

        for (int i=0; eIt!=eEndIt; ++eIt, i++) {

          GeometryType type = eIt->type();

          if (type.isHexahedron()) {

            const int hexaReordering[8] = {0, 1, 3, 2, 4, 5, 7, 6};
            for (int j=0; j<8; j++)
              dPtr[8*i + j] = indexSet.subIndex(*eIt, hexaReordering[j],dim)+1;

          } else if (type.isPrism()) {

            const int prismReordering[8] = {0, 1, 1, 2, 3, 4, 4, 5};
            for (int j=0; j<8; j++)
              dPtr[8*i + j] = indexSet.subIndex(*eIt, prismReordering[j],dim)+1;

          } else if (type.isPyramid()) {

            const int pyramidReordering[8] = {0, 1, 3, 2, 4, 4, 4, 4};
            for (int j=0; j<8; j++)
              dPtr[8*i + j] = indexSet.subIndex(*eIt, pyramidReordering[j], dim)+1;

          } else if (type.isTetrahedron()) {

            const int tetraReordering[8] = {0, 1, 2, 2, 3, 3, 3, 3};
            for (int j=0; j<8; j++)
              dPtr[8*i + j] = indexSet.subIndex(*eIt, tetraReordering[j],dim)+1;

          } else
            DUNE_THROW(NotImplemented, "Unknown element type encountered");

        }

      }

    } else {

      for (int i=0; eIt!=eEndIt; ++eIt, i++) {

        GeometryType type = eIt->type();

        if (type.isQuadrilateral()) {

          dPtr[i*4+0] = indexSet.subIndex(*eIt, 0, dim)+1;
          dPtr[i*4+1] = indexSet.subIndex(*eIt, 1, dim)+1;
          dPtr[i*4+2] = indexSet.subIndex(*eIt, 3, dim)+1;
          dPtr[i*4+3] = indexSet.subIndex(*eIt, 2, dim)+1;

        } else if (type.isTriangle()) {

          for (int j=0; j<3; j++)
            dPtr[i*maxVerticesPerElement+j] = indexSet.subIndex(*eIt, j, dim)+1;

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

  }

  // write material section to grid file

  AmiraMesh::Data* element_materials = new AmiraMesh::Data("Materials", elementLocation, McPrimType::mc_uint8, 1);
  amiramesh_.insert(element_materials);

  for(int i=0; i<noOfElements; i++)
    ((unsigned char*)element_materials->dataPtr())[i] = 0;

}

template<class GridView>
template<class GridType2>
void Dune::AmiraMeshWriter<GridView>::addLevelGrid(const GridType2& grid,
                                                   int level,
                                                   bool splitAll)
{
  addGrid(grid.levelView(level), splitAll);
}


template<class GridView>
template<class GridType2>
void Dune::AmiraMeshWriter<GridView>::addLeafGrid(const GridType2& grid, bool splitAll)
{
  addGrid(grid.leafView(), splitAll);
}


template<class GridView>
template<class DataContainer>
void Dune::AmiraMeshWriter<GridView>::addCellData(const DataContainer& data,
                                                  const GridView& gridView,
                                                  bool gridSplitUp)
{
  typedef typename GridView::template Codim<0>::Iterator ElementIterator;
  typedef typename GridView::Grid::ctype ct;
  typedef Dune::VirtualRefinement<dim, ct> Refinement;


  const typename GridView::IndexSet& indexSet = gridView.indexSet();

  //gridSplitUp tells programm that Amira "thinks" that all elements are tethraheda


  // Find out whether the grid contains only tetrahedra.  If yes, then
  // it is written in TetraGrid format.  If not, it is written in
  // hexagrid format (if gridSplitUp=false).
  bool containsOnlyTetrahedra =
    (indexSet.geomTypes(0).size()==1)
    && (indexSet.geomTypes(0)[0].isSimplex());

  // Get number of components
  const int ncomp = DataContainer::block_type::size;

  // Set the appropriate content type for 2D grid data, if no other
  // content type hasn't been set already
  if (dim==2 and amiramesh_.parameters.findBase("ContentType")==NULL)
    amiramesh_.parameters.set("ContentType", "HxTriangularData");

  if (!containsOnlyTetrahedra and dim==3 and !gridSplitUp)
  {
    AmiraMesh::Location* hexa_loc = new AmiraMesh::Location("Hexahedra", indexSet.size(0));
    amiramesh_.insert(hexa_loc);
  }

  AmiraMesh::Location* amLocation;

  if (data.size()==indexSet.size(0))
  {
    if (gridSplitUp || containsOnlyTetrahedra)
    {
      Dune::GeometryType coerceTo(Dune::GeometryType::simplex,dim);
      int noOfElements = 0;
      int count;

      for (int i=0; i<indexSet.geomTypes(0).size(); i++)
      {
        if (indexSet.geomTypes(0)[i].isSimplex())
          count=1;
        else
        {
          Refinement & refinement = Dune::buildRefinement<dim, ct>(indexSet.geomTypes(0)[i],coerceTo);
          count = refinement.nElements(0);
        }

        noOfElements += count * indexSet.size(indexSet.geomTypes(0)[i]);
      }

      amLocation = new AmiraMesh::Location((dim==2) ? "Triangles" : "Tetrahedra", noOfElements);
    }
    else
      amLocation = new AmiraMesh::Location("Hexahedra", data.size());

  }
  else
    DUNE_THROW(IOError, "AmiraMeshWriter::addCellData: BlockVector doesn't match the grid! For vertex based data use addVertexData");

  amiramesh_.insert(amLocation);

  // \todo Auto-detect data type
  AmiraMesh::Data* nodeData = new AmiraMesh::Data("Data", amLocation, McPrimType::mc_double, ncomp);
  amiramesh_.insert(nodeData);


  AmiraMesh::Field* nodeField;

  nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double, AmiraMesh::t_constant, nodeData);

  amiramesh_.insert(nodeField);

  // write the data into the AmiraMesh object
  typedef typename DataContainer::ConstIterator Iterator;

  Iterator dit    = data.begin();
  Iterator ditend = data.end();

  int i=0;
  if (gridSplitUp)
  {
    ElementIterator eIt    = gridView.template begin<0>();
    Dune::GeometryType coerceTo(Dune::GeometryType::simplex,dim);

    for (; dit!=ditend; ++dit)
    {
      Refinement & refinement = Dune::buildRefinement<dim, ct>(eIt->type(),coerceTo);
      int num_subsimplices=refinement.nElements(0);

      //Have to copy data if gridSplitUp because number_elements != number_data;
      for(int k=0; k<num_subsimplices; k++)
      {
        for (int j=0; j<ncomp; j++)
          ((double*)nodeData->dataPtr())[i++] = (*dit)[j];
      }
      ++eIt;
    }
  } else {

    // Write data directly
    for (; dit!=ditend; ++dit)
      for (int j=0; j<ncomp; j++)
        ((double*)nodeData->dataPtr())[i++] = (*dit)[j];

  }
}



template<class GridView>
template<class DataContainer>
void Dune::AmiraMeshWriter<GridView>::addVertexData(const DataContainer& data,
                                                    const GridView& gridView,
                                                    bool gridSplitUp)
{
  typedef typename GridView::template Codim<0>::Iterator ElementIterator;
  typedef typename GridView::Grid::ctype ct;
  typedef Dune::VirtualRefinement<dim, ct> Refinement;


  const typename GridView::IndexSet& indexSet = gridView.indexSet();

  //gridSplitUp tells programm that Amira "thinks" that all elements are tethraheda


  // Find out whether the grid contains only tetrahedra.  If yes, then
  // it is written in TetraGrid format.  If not, it is written in
  // hexagrid format (if gridSplitUp=false).
  bool containsOnlyTetrahedra =
    (indexSet.geomTypes(0).size()==1)
    && (indexSet.geomTypes(0)[0].isSimplex());

  // Get number of components
  const int ncomp = DataContainer::block_type::size;

  // Set the appropriate content type for 2D grid data, if no other
  // content type hasn't been set already
  if (dim==2 and amiramesh_.parameters.findBase("ContentType")==NULL)
    amiramesh_.parameters.set("ContentType", "HxTriangularData");

  if (!containsOnlyTetrahedra and dim==3 and !gridSplitUp)
  {
    AmiraMesh::Location* hexa_loc = new AmiraMesh::Location("Hexahedra", indexSet.size(0));
    amiramesh_.insert(hexa_loc);
  }

  AmiraMesh::Location* amLocation;

  if (data.size()==indexSet.size(dim))
  {
    // P1 data
    amLocation = new AmiraMesh::Location("Nodes", data.size());
  }
  else
    DUNE_THROW(IOError, "AmiraMeshWriter::addVertexData: BlockVector doesn't match the grid! For element based data use addCellData");

  amiramesh_.insert(amLocation);

  // \todo Auto-detect data type
  AmiraMesh::Data* nodeData = new AmiraMesh::Data("Data", amLocation, McPrimType::mc_double, ncomp);
  amiramesh_.insert(nodeData);


  AmiraMesh::Field* nodeField;

  if (containsOnlyTetrahedra || dim==2 || gridSplitUp)
  {
    nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double, AmiraMesh::t_linear, nodeData);
  }
  else
  {
    nodeField = new AmiraMesh::Field("sol", ncomp, McPrimType::mc_double, AmiraMesh::t_trilinear, nodeData);
  }

  amiramesh_.insert(nodeField);

  // write the data into the AmiraMesh object
  typedef typename DataContainer::ConstIterator Iterator;

  Iterator dit    = data.begin();
  Iterator ditend = data.end();

  int i=0;

  for (; dit!=ditend; ++dit)
  {
    for (int j=0; j<ncomp; j++)
      ((double*)nodeData->dataPtr())[i++] = (*dit)[j];
  }

}


template<class GridView>
void Dune::AmiraMeshWriter<GridView>::write(const std::string& filename,
                                            bool ascii) const
{
  // Actually write the file
  if(!amiramesh_.write(filename.c_str(), ascii))
    DUNE_THROW(IOError, "Writing geometry file '" << filename << "' failed!");

  std::cout << "Grid written successfully to: " << filename << std::endl;
}


template<class GridView>
template<class DataContainer>
void Dune::AmiraMeshWriter<GridView>::addUniformData(const GridView& gridView,
                                                     const array<unsigned int, dim>& n,
                                                     const DataContainer& data)
{
  dune_static_assert(dim==2 || dim==3, "You can only write 2d and 3d uniform data to AmiraMesh");

  // ///////////////////////////////////////////
  //   Detect grid bounding box
  // ///////////////////////////////////////////
  float bbox[2*dim];
  for (int i=0; i<dim; i++) {
    bbox[2*i  ] =  std::numeric_limits<double>::max();
    bbox[2*i+1] = -std::numeric_limits<double>::max();
  }

  typename GridView::template Codim<dim>::Iterator vIt    = gridView.template begin<dim>();
  typename GridView::template Codim<dim>::Iterator vEndIt = gridView.template end<dim>();

  for (; vIt!=vEndIt; ++vIt)
    for (int i=0; i<dim; i++) {
      bbox[2*i]   = std::min((double)bbox[2*i],   vIt->geometry().corner(0)[i]);
      bbox[2*i+1] = std::max((double)bbox[2*i+1], vIt->geometry().corner(0)[i]);
    }

  // Set the appropriate content type
  if (dim==2)
    amiramesh_.parameters.set("ContentType", "HxField2d");

  amiramesh_.parameters.set("BoundingBox", 2*dim, bbox);
  amiramesh_.parameters.set("CoordType", "uniform");

  AmiraMesh::Location* loc = amiramesh_.findLocation("Lattice");
  int dims[dim];
  for (int i=0; i<dim; i++)
    dims[i] = n[i];

  if (!loc) {
    loc = new AmiraMesh::Location("Lattice", dim, dims);
    amiramesh_.insert(loc);
  }

  // set up data
  // we assume that at least the inner vector follows ISTL conventions
  // unfortunately istl and std containers are not compatible
  // in that it is not possible to extract the size
  // of array and FieldVector by the same method
  int numComponents = DataContainer::value_type::size;

  AmiraMesh::Data* amData =
    new AmiraMesh::Data("Data", loc, McPrimType::mc_double, numComponents);
  amiramesh_.insert(amData);

  // ////////////////////////////////////////////////////////
  //   Write the data
  // ////////////////////////////////////////////////////////
  typedef typename DataContainer::const_iterator iterator;
  iterator it    = data.begin();
  iterator endIt = data.end();

  int i=0;
  for (; it!=endIt; ++it)
    for (int j=0; j<numComponents; j++, i++)
      ((double*)amData->dataPtr())[i] = (*it)[j];

}
