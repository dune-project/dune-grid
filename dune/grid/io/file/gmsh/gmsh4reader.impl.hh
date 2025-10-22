// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <array>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <set>

#include <dune/grid/io/file/gmsh/utility/errors.hh>
#include <dune/grid/io/file/gmsh/utility/string.hh>

namespace Dune::Impl::Gmsh
{

  template <class G, class C, class S>
  std::map<int, S> Gmsh4Reader<G,C,S>::elementType_ {
    {1, 2u},  // 2-node line.
    {2, 3u},  // 3-node triangle.
    {3, 4u},  // 4-node quadrangle.
    {4, 4u},  // 4-node tetrahedron.
    {5, 8u},  // 8-node hexahedron.
    {6, 6u},  // 6-node prism.
    {7, 5u},  // 5-node pyramid.
    {8, 3u},  // 3-node second order line.
    {9, 6u},  // 6-node second order triangle.
    {10, 9u}, // 9-node second order quadrangle.
    {11, 10u}, // 10-node second order tetrahedron.
    {12, 27u}, // 27-node second order hexahedron.
    {13, 18u}, // 18-node second order prism.
    {14, 14u}, // 14-node second order pyramid.
    {15, 1u}, // 1-node point.
    {16, 8u}, // 8-node second order quadrangle.
    {17, 20u}, // 20-node second order hexahedron.
    {18, 15u}, // 15-node second order prism.
    {19, 13u}, // 13-node second order pyramid.
    {20, 9u}, // 9-node third order incomplete triangle.
    {21, 10u}, // 10-node third order triangle.
    {22, 12u}, // 12-node fourth order incomplete triangle.
    {23, 15u}, // 15-node fourth order triangle.
    {24, 15u}, // 15-node fifth order incomplete triangle.
    {25, 21u}, // 21-node fifth order complete triangle.
    {26, 4u}, // 4-node third order line.
    {27, 5u}, // 5-node fourth order line.
    {28, 6u}, // 6-node fifth order line.
    {29, 20u}, // 20-node third order tetrahedron.
    {30, 35u}, // 35-node fourth order tetrahedron.
    {31, 56u}, // 56-node fifth order tetrahedron.
    {32, 22u}, // 22-node tetrahedron.
    {33, 28u}, // 28-node tetrahedron.
    //{34, ?},  // polygon.
    //{35, ?},  // polygon.
    {36, 16u}, // 16-node quadrangle.
    {37, 25u}, // 25-node quadrangle.
    {38, 36u}, // 36-node quadrangle.
    {39, 12u}, // 12-node quadrangle.
    {40, 16u}, // 16-node quadrangle (I).
    {41, 20u}, // 20-node quadrangle.
    {42, 28u}, // 28-node triangle.
    {43, 36u}, // 36-node triangle.
    {44, 45u}, // 45-node triangle.
    {45, 55u}, // 55-node triangle.
    {46, 66u}, // 66-node triangle.
    {47, 49u}, // 49-node quadrangle.
    {48, 64u}, // 64-node quadrangle.
    {49, 81u}, // 81-node quadrangle.
    {50, 100u}, // 100-node quadrangle.
    {51, 121u}, // 121-node quadrangle.
    {52, 18u}, // 18-node triangle.
    {53, 21u}, // 21-node triangle (I).
    {54, 24u}, // 24-node triangle.
    {55, 27u}, // 27-node triangle.
    {56, 30u}, // 30-node triangle.
    {57, 24u}, // 24-node quadrangle.
    {58, 28u}, // 28-node quadrangle.
    {59, 32u}, // 32-node quadrangle.
    {60, 36u}, // 36-node quadrangle (I).
    {61, 40u}, // 40-node quadrangle.
    {62, 7u}, // 7-node line.
    {63, 8u}, // 8-node line.
    {64, 9u}, // 9-node line.
    {65, 10u}, // 10-node line.
    {66, 11u}, // 11-node line.
    //{67, ?},  // line.
    //{68, ?},  // triangle.
    //{69, ?},  // polygon.
    //{70, ?},  // line.
    {71, 84u}, // 84-node tetrahedron.
    {72, 120u}, // 120-node tetrahedron.
    {73, 165u}, // 165-node tetrahedron.
    {74, 220u}, // 220-node tetrahedron.
    {75, 286u}, // 286-node tetrahedron.
    {79, 34u}, // 34-node incomplete tetrahedron.
    {80, 40u}, // 40-node incomplete tetrahedron.
    {81, 46u}, // 46-node incomplete tetrahedron.
    {82, 52u}, // 52-node incomplete tetrahedron.
    {83, 58u}, // 58-node incomplete tetrahedron.
    {84, 1u}, // 1-node line.
    {85, 1u}, // 1-node triangle.
    {86, 1u}, // 1-node quadrangle.
    {87, 1u}, // 1-node tetrahedron.
    {88, 1u}, // 1-node hexahedron.
    {89, 1u}, // 1-node prism.
    {90, 40u}, // 40-node prism.
    {91, 75u}, // 75-node prism.
    {92, 64u}, // 64-node third order hexahedron.
    {93, 125u}, // 125-node fourth order hexahedron.
    {94, 216u}, // 216-node hexahedron.
    {95, 343u}, // 343-node hexahedron.
    {96, 512u}, // 512-node hexahedron.
    {97, 729u}, // 729-node hexahedron.
    {98, 1000u},// 1000-node hexahedron.
    {99, 32u}, // 32-node incomplete hexahedron.
    {100, 44u}, // 44-node incomplete hexahedron.
    {101, 56u}, // 56-node incomplete hexahedron.
    {102, 68u}, // 68-node incomplete hexahedron.
    {103, 80u}, // 80-node incomplete hexahedron.
    {104, 92u}, // 92-node incomplete hexahedron.
    {105, 104u},// 104-node incomplete hexahedron.
    {106, 126u},// 126-node prism.
    {107, 196u},// 196-node prism.
    {108, 288u},// 288-node prism.
    {109, 405u},// 405-node prism.
    {110, 550u},// 550-node prism.
    {111, 24u}, // 24-node incomplete prism.
    {112, 33u}, // 33-node incomplete prism.
    {113, 42u}, // 42-node incomplete prism.
    {114, 51u}, // 51-node incomplete prism.
    {115, 60u}, // 60-node incomplete prism.
    {116, 69u}, // 69-node incomplete prism.
    {117, 78u}, // 78-node incomplete prism.
    {118, 30u}, // 30-node pyramid.
    {119, 55u}, // 55-node pyramid.
    {120, 91u}, // 91-node pyramid.
    {121, 140u},// 140-node pyramid.
    {122, 204u},// 204-node pyramid.
    {123, 285u},// 285-node pyramid.
    {124, 385u},// 385-node pyramid.
    {125, 21u}, // 21-node incomplete pyramid.
    {126, 29u}, // 29-node incomplete pyramid.
    {127, 37u}, // 37-node incomplete pyramid.
    {128, 45u}, // 45-node incomplete pyramid.
    {129, 53u}, // 53-node incomplete pyramid.
    {130, 61u}, // 61-node incomplete pyramid.
    {131, 69u}, // 69-node incomplete pyramid.
    {132, 1u}, // 1-node pyramid.
    //{133, ?}, // point.
    //{134, ?}, // line.
    //{135, ?}, // triangle.
    //{136, ?}, // tetrahedron.
    {137, 16u}, // 16-node tetrahedron.
    //{138, ?}, // triangle (mini).
    //{139, ?}, // tetrahedron (mini).
    {140, 4u}, // 4-node triangle.
  };

  template <class G, class C, class S>
  std::map<std::string, typename Gmsh4Reader<G,C,S>::Sections> Gmsh4Reader<G,C,S>::sections_ {
    {"MeshFormat", Sections::MESH_FORMAT},
    {"PhysicalNames", Sections::PHYSICAL_NAMES},
    {"Entities", Sections::ENTITIES},
    {"PartitionedEntities", Sections::PARTITIONED_ENTITIES},
    {"Nodes", Sections::NODES},
    {"Elements", Sections::ELEMENTS},
    {"Periodic", Sections::PERIODIC},
    {"GhostElements", Sections::GHOST_ELEMENTS},
    {"Parametrization", Sections::PARAMETRIZATION},
    {"NodeData", Sections::NODE_DATA},
    {"ElementData", Sections::ELEMENT_DATA},
    {"ElementNodeData", Sections::ELEMENT_NODE_DATA},
    {"InterpolationScheme", Sections::INTERPOLATION_SCHEME}
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::PhysicalNamesAttributes
  {
    int dim;
    int tag;
    std::string name;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::PointAttributes
  {
    int tag;
    std::array<double,3> xyz;
    std::vector<int> physicalTags;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::EntityAttributes
  {
    int tag;
    std::array<double,3> min_xyz;
    std::array<double,3> max_xyz;
    std::vector<int> physicalTags;
    std::vector<int> boundingEntities;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::GhostAttributes
  {
    int tag;
    int partition;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::NodeAttributes
  {
    struct Node
    {
      size_type tag;
      std::array<double,3> xyz;
      std::array<double,3> uvw;
    };

    int entityDim;
    int entityTag;
    int parametric;
    std::vector<Node> nodes;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::ElementAttributes
  {
    struct Element
    {
      size_type tag;
      std::vector<size_type> nodes;
    };

    int entityDim;
    int entityTag;
    int elementType;
    std::vector<Element> elements;
  };

  template <class G, class C, class S>
  struct Gmsh4Reader<G,C,S>::PeriodicAttributes
  {
    struct Association
    {
      size_type tag;
      size_type tagMaster;
    };

    int entityDim;
    int entityTag;
    int entityTagMaster;
    std::vector<double> affine;
    std::vector<Association> correspondingNodes;
  };


  template <class G, class C, class S>
  template <class T>
  void Gmsh4Reader<G,C,S>::readValueBinary(std::ifstream& input, T &v)
  {
    const std::size_t size = sizeof(T);
    input.read(reinterpret_cast<char*>(&v), size);
    if (swap)
    {
      char* vBinary = reinterpret_cast<char*>(&v);
      for (std::size_t i=0; i<size/2; ++i)
        std::swap(vBinary[i],vBinary[size-1-i]);
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::read (std::string const& filename, bool fillCreator)
  {
    const std::filesystem::path filePath(filename);

    // check whether file exists!
    if (!std::filesystem::exists(filePath))
      DUNE_THROW(IOError, "File " << filename << " does not exist!");

    std::ifstream input(filename, std::ios_base::in | std::ios_base::binary);
    GMSH4_ASSERT(input.is_open());

    const std::string ext = filePath.extension().string();
    if (ext == ".msh") {

      // Get first line
      std::string line;
      std::getline(input, line);
      ltrim(line);

      // Read the mesh format section -- it must make up the first two lines
      GMSH4_ASSERT_MSG((line=="$MeshFormat"), "First line of file is not $MeshFormat!");

      double version = 0.0;
      int file_type = 0;
      int data_size = 0;
      readMeshFormat(input, version, file_type, data_size);

      // Test whether version is supported
      GMSH4_ASSERT_MSG((version >= 4.0 && version < 5.0),
                       "Can only read gmsh files versions >= 4.0 and < 5.0");

      // Further sanity checking
      GMSH4_ASSERT_MSG(file_type == 0 || file_type == 1,
                       "Invalid file-type: 0 for ASCII mode, 1 for binary mode");
      GMSH4_ASSERT_MSG(data_size >= 4 && data_size <= 16,
                       "Invalid data-size range: should be in {4, 16}");
      GMSH4_ASSERT_MSG(file_type != 1 || data_size == sizeof(size_type),
                       "Invalid data-size: must be sizeof(size_t)");

      // Rewind the stream, so readSerialFileFromStream can start from the top again
      input.clear();
      input.seekg(0, std::ios::beg);

      readSerialFileFromStream(input, fillCreator);
      pieces_.push_back(filename);
    } else if (ext == ".pro") {
      readParallelFileFromStream(input, comm().rank(), comm().size(), fillCreator);
    } else {
      DUNE_THROW(IOError, "File has unknown file-extension '" << ext << "'. Allowed are only '.msh' and '.pro'.");
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readSerialFileFromStream (std::ifstream& input, bool fillCreator)
  {
    clear();

    // MeshFormat section
    int file_type = 0;
    int data_size = 0;

    Sections section = Sections::NO_SECTION;
    for (std::string line; std::getline(input, line); ) {
      ltrim(line);

      // detect current section
      for (auto const& s : sections_) {
        if (isSection(line, s.first, section)) {
          section = s.second;
          break;
        } else if (isSection(line, "End" + s.first, section, s.second)) {
          section = Sections::NO_SECTION;
          break;
        }
      }

      switch (section) {
      case Sections::MESH_FORMAT: {
        double version = 0.0;
        readMeshFormat(input, version, file_type, data_size);
        break;
      }
      case Sections::PHYSICAL_NAMES:
        readPhysicalNames(input); break;
      case Sections::ENTITIES:
        if(file_type == 0) readEntitiesAscii(input);
        else readEntitiesBinary(input);
        break;
      case Sections::PARTITIONED_ENTITIES:
        if(file_type == 0) readPartitionedEntitiesAscii(input);
        else readPartitionedEntitiesBinary(input);
        break;
      case Sections::NODES:
        if(file_type == 0) readNodesAscii(input);
        else readNodesBinary(input);
        break;
      case Sections::ELEMENTS:
        if(file_type == 0) readElementsAscii(input);
        else readElementsBinary(input);
        break;
      case Sections::PERIODIC:
        readPeriodic(input); break;
      case Sections::GHOST_ELEMENTS:
        readGhostElements(input); break;
      case Sections::PARAMETRIZATION:
        readParametrization(input); break;
      case Sections::NODE_DATA:
        readNodeData(input); break;
      case Sections::ELEMENT_DATA:
        readElementData(input); break;
      case Sections::ELEMENT_NODE_DATA:
        readElementNodeData(input); break;
      case Sections::INTERPOLATION_SCHEME:
        readInterpolationScheme(input); break;
      default:
        // do nothing
        break;
      }
    }

    if (fillCreator)
      fillGridCreator();
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readParallelFileFromStream (std::ifstream& input, int commRank, int commSize, bool fillCreator)
  {
    clear();
    DUNE_THROW(Dune::NotImplemented, "Reading parallel .pro files not yet implemented.");
    if (fillCreator)
      fillGridCreator();
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readMeshFormat (std::ifstream& input, double& version, int& file_type, int& data_size)
  {
    std::string line;
    std::getline(input, line);
    std::istringstream stream(line);
    stream >> version >> file_type >> data_size;
    if (file_type != 0) {
      int one;
      input.read(reinterpret_cast<char*>(&one), sizeof(int));
      if (one != 1) swap = true;
      std::getline(input, line);
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readPhysicalNames (std::ifstream& input)
  {
    int numPhysicalNames = 0;
    std::string line;
    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numPhysicalNames;
    }

    for (int i = 0; i < numPhysicalNames; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PhysicalNamesAttributes attr;
      stream >> attr.dim >> attr.tag;
      readString(stream, attr.name);

      physicalNames_.push_back(attr);
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readEntitiesAscii (std::ifstream& input)
  {
    size_type numPoints = 0, numCurves = 0, numSurfaces = 0, numVolumes = 0;

    std::string line;
    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numPoints >> numCurves >> numSurfaces >> numVolumes;
    }

    // points
    points_.reserve(numPoints);
    for (size_type i = 0; i < numPoints; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PointAttributes attr;
      stream >> attr.tag >> attr.xyz[0] >> attr.xyz[1] >> attr.xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;

      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      points_.push_back(attr);
    }
    GMSH4_ASSERT(points_.size() == numPoints);

    // curves
    curves_.reserve(numCurves);
    for (size_type i = 0; i < numCurves; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      EntityAttributes attr;
      stream >> attr.tag
      >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingPoints = 0;
      stream >> numBoundingPoints;
      attr.boundingEntities.resize(numBoundingPoints);
      for (size_type j = 0; j < numBoundingPoints; ++j)
        stream >> attr.boundingEntities[j];

      curves_.push_back(attr);
    }
    GMSH4_ASSERT(curves_.size() == numCurves);

    // surfaces
    surfaces_.reserve(numSurfaces);
    for (size_type i = 0; i < numSurfaces; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      EntityAttributes attr;
      stream >> attr.tag
      >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingCurves = 0;
      stream >> numBoundingCurves;
      attr.boundingEntities.resize(numBoundingCurves);
      for (size_type j = 0; j < numBoundingCurves; ++j)
        stream >> attr.boundingEntities[j];

      surfaces_.push_back(attr);
    }
    GMSH4_ASSERT(surfaces_.size() == numSurfaces);

    // volumes
    volumes_.reserve(numVolumes);
    for (size_type i = 0; i < numVolumes; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      EntityAttributes attr;
      stream >> attr.tag
      >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingSurfaces = 0;
      stream >> numBoundingSurfaces;
      attr.boundingEntities.resize(numBoundingSurfaces);
      for (size_type j = 0; j < numBoundingSurfaces; ++j)
        stream >> attr.boundingEntities[j];

      volumes_.push_back(attr);
    }
    GMSH4_ASSERT(volumes_.size() == numVolumes);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readEntitiesBinary (std::ifstream& input)
  {
    size_type numPoints = 0, numCurves = 0, numSurfaces = 0, numVolumes = 0;

    readValueBinary(input, numPoints);
    readValueBinary(input, numCurves);
    readValueBinary(input, numSurfaces);
    readValueBinary(input, numVolumes);

    // points
    points_.reserve(numPoints);
    for (size_type i = 0; i < numPoints; ++i) {
      PointAttributes attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.xyz[0]);
      readValueBinary(input, attr.xyz[1]);
      readValueBinary(input, attr.xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      points_.push_back(attr);
    }
    GMSH4_ASSERT(points_.size() == numPoints);

    // curves
    curves_.reserve(numCurves);
    for (size_type i = 0; i < numCurves; ++i) {
      EntityAttributes attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingPoints = 0;
      readValueBinary(input, numBoundingPoints);
      attr.boundingEntities.resize(numBoundingPoints);
      for (size_type j = 0; j < numBoundingPoints; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      curves_.push_back(attr);
    }
    GMSH4_ASSERT(curves_.size() == numCurves);

    // surfaces
    surfaces_.reserve(numSurfaces);
    for (size_type i = 0; i < numSurfaces; ++i) {
      EntityAttributes attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingCurves = 0;
      readValueBinary(input, numBoundingCurves);
      attr.boundingEntities.resize(numBoundingCurves);
      for (size_type j = 0; j < numBoundingCurves; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      surfaces_.push_back(attr);
    }
    GMSH4_ASSERT(surfaces_.size() == numSurfaces);

    // volumes
    volumes_.reserve(numVolumes);
    for (size_type i = 0; i < numVolumes; ++i) {
      EntityAttributes attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingSurfaces = 0;
      readValueBinary(input, numBoundingSurfaces);
      attr.boundingEntities.resize(numBoundingSurfaces);
      for (size_type j = 0; j < numBoundingSurfaces; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      volumes_.push_back(attr);
    }
    GMSH4_ASSERT(volumes_.size() == numVolumes);
    std::string line;
    std::getline(input, line);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readPartitionedEntitiesAscii (std::ifstream& input)
  {
    size_type numGhostEntities = 0;
    size_type numPoints = 0, numCurves = 0, numSurfaces = 0, numVolumes = 0;

    std::string line;
    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numPartitions_;
    }

    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numGhostEntities;

      // ghost entities
      ghostEntities_.reserve(numGhostEntities);
      for (size_type i = 0; i < numGhostEntities; ++i) {
        GhostAttributes attr;
        stream >> attr.tag >> attr.partition;

        ghostEntities_.push_back(attr);
      }
      GMSH4_ASSERT(ghostEntities_.size() == numGhostEntities);
    }

    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numPoints >> numCurves >> numSurfaces >> numVolumes;
    }

    // points
    partitionedPoints_.reserve(numPoints);
    for (size_type i = 0; i < numPoints; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PartitionedAttributes<PointAttributes> attr;
      stream >> attr.tag >> attr.parentDim >> attr.parentTag;

      size_type numPartitions = 0;
      stream >> numPartitions;
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        stream >> attr.partitions[j];

      stream >> attr.xyz[0] >> attr.xyz[1] >> attr.xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      partitionedPoints_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedPoints_.size() == numPoints);

    // curves
    partitionedCurves_.reserve(numCurves);
    for (size_type i = 0; i < numCurves; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PartitionedAttributes<EntityAttributes> attr;
      stream >> attr.tag >> attr.parentDim >> attr.parentTag;

      size_type numPartitions = 0;
      stream >> numPartitions;
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        stream >> attr.partitions[j];

      stream >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingPoints = 0;
      stream >> numBoundingPoints;
      attr.boundingEntities.resize(numBoundingPoints);
      for (size_type j = 0; j < numBoundingPoints; ++j)
        stream >> attr.boundingEntities[j];

      partitionedCurves_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedCurves_.size() == numCurves);

    // surfaces
    partitionedSurfaces_.reserve(numSurfaces);
    for (size_type i = 0; i < numSurfaces; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PartitionedAttributes<EntityAttributes> attr;
      stream >> attr.tag >> attr.parentDim >> attr.parentTag;

      size_type numPartitions = 0;
      stream >> numPartitions;
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        stream >> attr.partitions[j];

      stream >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingCurves = 0;
      stream >> numBoundingCurves;
      attr.boundingEntities.resize(numBoundingCurves);
      for (size_type j = 0; j < numBoundingCurves; ++j)
        stream >> attr.boundingEntities[j];

      partitionedSurfaces_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedSurfaces_.size() == numSurfaces);

    // volumes
    partitionedVolumes_.reserve(numVolumes);
    for (size_type i = 0; i < numVolumes; ++i) {
      if (!std::getline(input,line))
        break;

      std::istringstream stream(line);
      PartitionedAttributes<EntityAttributes> attr;
      stream >> attr.tag >> attr.parentDim >> attr.parentTag;

      size_type numPartitions = 0;
      stream >> numPartitions;
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        stream >> attr.partitions[j];

      stream >> attr.min_xyz[0] >> attr.min_xyz[1] >> attr.min_xyz[2]
      >> attr.max_xyz[0] >> attr.max_xyz[1] >> attr.max_xyz[2];

      size_type numPhysicalTags = 0;
      stream >> numPhysicalTags;
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        stream >> attr.physicalTags[j];

      size_type numBoundingSurfaces = 0;
      stream >> numBoundingSurfaces;
      attr.boundingEntities.resize(numBoundingSurfaces);
      for (size_type j = 0; j < numBoundingSurfaces; ++j)
        stream >> attr.boundingEntities[j];

      partitionedVolumes_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedVolumes_.size() == numVolumes);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readPartitionedEntitiesBinary (std::ifstream& input)
  {
    size_type numGhostEntities = 0;
    size_type numPoints = 0, numCurves = 0, numSurfaces = 0, numVolumes = 0;

    readValueBinary(input, numPartitions_);

    // ghost entities
    readValueBinary(input, numGhostEntities);
    ghostEntities_.reserve(numGhostEntities);
    for (size_type i = 0; i < numGhostEntities; ++i) {
      GhostAttributes attr;
      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.partition);
      ghostEntities_.push_back(attr);
    }
    GMSH4_ASSERT(ghostEntities_.size() == numGhostEntities);

    readValueBinary(input, numPoints);
    readValueBinary(input, numCurves);
    readValueBinary(input, numSurfaces);
    readValueBinary(input, numVolumes);

    // points
    partitionedPoints_.reserve(numPoints);
    for (size_type i = 0; i < numPoints; ++i) {
      PartitionedAttributes<PointAttributes> attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.parentDim);
      readValueBinary(input, attr.parentTag);

      size_type numPartitions = 0;
      readValueBinary(input, numPartitions);
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        readValueBinary(input, attr.partitions[j]);

      readValueBinary(input, attr.xyz[0]);
      readValueBinary(input, attr.xyz[1]);
      readValueBinary(input, attr.xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      partitionedPoints_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedPoints_.size() == numPoints);

    // curves
    partitionedCurves_.reserve(numCurves);
    for (size_type i = 0; i < numCurves; ++i) {
      PartitionedAttributes<EntityAttributes> attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.parentDim);
      readValueBinary(input, attr.parentTag);

      size_type numPartitions = 0;
      readValueBinary(input, numPartitions);
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        readValueBinary(input, attr.partitions[j]);

      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingPoints = 0;
      readValueBinary(input, numBoundingPoints);
      attr.boundingEntities.resize(numBoundingPoints);
      for (size_type j = 0; j < numBoundingPoints; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      partitionedCurves_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedCurves_.size() == numCurves);

    // surfaces
    partitionedSurfaces_.reserve(numSurfaces);
    for (size_type i = 0; i < numSurfaces; ++i) {
      PartitionedAttributes<EntityAttributes> attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.parentDim);
      readValueBinary(input, attr.parentTag);

      size_type numPartitions = 0;
      readValueBinary(input, numPartitions);
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        readValueBinary(input, attr.partitions[j]);

      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingCurves = 0;
      readValueBinary(input, numBoundingCurves);
      attr.boundingEntities.resize(numBoundingCurves);
      for (size_type j = 0; j < numBoundingCurves; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      partitionedSurfaces_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedSurfaces_.size() == numSurfaces);

    // volumes
    partitionedVolumes_.reserve(numVolumes);
    for (size_type i = 0; i < numVolumes; ++i) {
      PartitionedAttributes<EntityAttributes> attr;

      readValueBinary(input, attr.tag);
      readValueBinary(input, attr.parentDim);
      readValueBinary(input, attr.parentTag);

      size_type numPartitions = 0;
      readValueBinary(input, numPartitions);
      attr.partitions.resize(numPartitions);
      for (size_type j = 0; j < numPartitions; ++j)
        readValueBinary(input, attr.partitions[j]);

      readValueBinary(input, attr.min_xyz[0]);
      readValueBinary(input, attr.min_xyz[1]);
      readValueBinary(input, attr.min_xyz[2]);
      readValueBinary(input, attr.max_xyz[0]);
      readValueBinary(input, attr.max_xyz[1]);
      readValueBinary(input, attr.max_xyz[2]);

      size_type numPhysicalTags = 0;
      readValueBinary(input, numPhysicalTags);
      attr.physicalTags.resize(numPhysicalTags);
      for (size_type j = 0; j < numPhysicalTags; ++j)
        readValueBinary(input, attr.physicalTags[j]);

      size_type numBoundingSurfaces = 0;
      readValueBinary(input, numBoundingSurfaces);
      attr.boundingEntities.resize(numBoundingSurfaces);
      for (size_type j = 0; j < numBoundingSurfaces; ++j)
        readValueBinary(input, attr.boundingEntities[j]);

      partitionedVolumes_.push_back(attr);
    }
    GMSH4_ASSERT(partitionedVolumes_.size() == numVolumes);
    std::string line;
    std::getline(input, line);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readNodesAscii (std::ifstream& input)
  {
    size_type numEntityBlocks = 0;

    std::string line;
    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numEntityBlocks >> numNodes_ >> minNodeTag_ >> maxNodeTag_;
    }

    nodes_.resize(numEntityBlocks);
    for (size_type i = 0; i < numEntityBlocks; ++i) {
      if (!std::getline(input,line))
        break;

      auto& entityBlock = nodes_[i];
      size_type numNodesInBlock = 0;

      std::istringstream stream(line);
      stream >> entityBlock.entityDim >> entityBlock.entityTag >> entityBlock.parametric >> numNodesInBlock;

      entityBlock.nodes.resize(numNodesInBlock);
      for (size_type j = 0; j < numNodesInBlock; ++j) {
        if (!std::getline(input,line))
          break;

        char* end;
        entityBlock.nodes[j].tag = std::strtoul(line.data(), &end, 10);
      }

      for (size_type j = 0; j < numNodesInBlock; ++j) {
        if (!std::getline(input,line))
          break;

        auto& node = entityBlock.nodes[j];

        std::istringstream stream(line);
        stream >> node.xyz[0] >> node.xyz[1] >> node.xyz[2];

        if (entityBlock.parametric && entityBlock.entityDim >= 1)
          stream >> node.uvw[0];
        if (entityBlock.parametric && entityBlock.entityDim >= 2)
          stream >> node.uvw[1];
        if (entityBlock.parametric && entityBlock.entityDim == 3)
          stream >> node.uvw[2];
      }
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readNodesBinary (std::ifstream& input)
  {
    size_type numEntityBlocks = 0;
    readValueBinary(input, numEntityBlocks);
    readValueBinary(input, numNodes_);
    readValueBinary(input, minNodeTag_);
    readValueBinary(input, maxNodeTag_);

    nodes_.resize(numEntityBlocks);
    for (size_type i = 0; i < numEntityBlocks; ++i) {
      auto& entityBlock = nodes_[i];
      size_type numNodesInBlock = 0;

      readValueBinary(input, entityBlock.entityDim);
      readValueBinary(input, entityBlock.entityTag);
      readValueBinary(input, entityBlock.parametric);
      readValueBinary(input, numNodesInBlock);

      entityBlock.nodes.resize(numNodesInBlock);
      for (size_type j = 0; j < numNodesInBlock; ++j)
        readValueBinary(input, entityBlock.nodes[j].tag);

      for (size_type j = 0; j < numNodesInBlock; ++j) {
        auto& node = entityBlock.nodes[j];

        readValueBinary(input, node.xyz[0]);
        readValueBinary(input, node.xyz[1]);
        readValueBinary(input, node.xyz[2]);

        if (entityBlock.parametric && entityBlock.entityDim >= 1)
          readValueBinary(input, node.uvw[0]);
        if (entityBlock.parametric && entityBlock.entityDim >= 2)
          readValueBinary(input, node.uvw[1]);
        if (entityBlock.parametric && entityBlock.entityDim == 3)
          readValueBinary(input, node.uvw[2]);
      }
    }
    std::string line;
    std::getline(input, line);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readElementsAscii (std::ifstream& input)
  {
    size_type numEntityBlocks = 0;
    size_type numElementsInBlock = 0;

    std::string line;
    std::getline(input, line);
    {
      std::istringstream stream(line);
      stream >> numEntityBlocks >> numElements_ >> minElementTag_ >> maxElementTag_;
    }

    elements_.resize(numEntityBlocks);
    for (size_type i = 0; i < numEntityBlocks; ++i) {
      if (!std::getline(input,line))
        break;

      auto& entityBlock = elements_[i];

      std::istringstream stream(line);
      stream >> entityBlock.entityDim >> entityBlock.entityTag >> entityBlock.elementType >> numElementsInBlock;

      size_type numNodes = elementType_[entityBlock.elementType];

      entityBlock.elements.resize(numElementsInBlock);
      for (size_type j = 0; j < numElementsInBlock; ++j) {
        if (!std::getline(input,line))
          break;

        auto& element = entityBlock.elements[j];

        std::istringstream stream(line);
        stream >> element.tag;
        element.nodes.resize(numNodes);
        for (size_type k = 0; k < numNodes; ++k)
          stream >> element.nodes[k];
      }
    }
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readElementsBinary (std::ifstream& input)
  {
    size_type numEntityBlocks = 0;
    size_type numElementsInBlock = 0;

    readValueBinary(input, numEntityBlocks);
    readValueBinary(input, numElements_);
    readValueBinary(input, minElementTag_);
    readValueBinary(input, maxElementTag_);

    elements_.resize(numEntityBlocks);
    for (size_type i = 0; i < numEntityBlocks; ++i) {
      auto& entityBlock = elements_[i];

      readValueBinary(input, entityBlock.entityDim);
      readValueBinary(input, entityBlock.entityTag);
      readValueBinary(input, entityBlock.elementType);
      readValueBinary(input, numElementsInBlock);

      size_type numNodes = elementType_[entityBlock.elementType];

      entityBlock.elements.resize(numElementsInBlock);
      for (size_type j = 0; j < numElementsInBlock; ++j) {
        auto& element = entityBlock.elements[j];
        readValueBinary(input, element.tag);
        element.nodes.resize(numNodes);
        for (size_type k = 0; k < numNodes; ++k)
          readValueBinary(input, element.nodes[k]);
      }
    }
    std::string line;
    std::getline(input, line);
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readPeriodic (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readPeriodic() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readGhostElements (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readGhostElements() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readParametrization (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readParametrization() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readNodeData (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readNodeData() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readElementData (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readElementData() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readElementNodeData (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readElementNodeData() is not yet implemented. Section will be ignored." << std::endl;
  }

  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::readInterpolationScheme (std::ifstream& /*input*/)
  {
    std::cout << "WARNING: readInterpolationScheme() is not yet implemented. Section will be ignored." << std::endl;
  }


  template <class G, class C, class S>
  void Gmsh4Reader<G,C,S>::fillGridCreator (bool insertPieces)
  {
    if (!nodes_.empty())
      creator_->insertVertices(numNodes_, {minNodeTag_, maxNodeTag_}, nodes_);
    if (!elements_.empty()) {
      std::set<int> boundaryEntities;
      creator_->insertElements(numElements_, {minElementTag_, maxElementTag_}, elements_, boundaryEntities);
    }
    if (insertPieces)
      creator_->insertPieces(pieces_);
  }

} // end namespace Dune::Impl::Gmsh
