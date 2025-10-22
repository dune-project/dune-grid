// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_GMSH4READER_HH
#define DUNE_GRID_IO_FILE_GMSH_GMSH4READER_HH

#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/gmsh/filereader.hh>
#include <dune/grid/io/file/gmsh/gridcreators/continuousgridcreator.hh> // default GridCreator

namespace Dune::Impl::Gmsh
{
  /// File-Reader for GMsh unstructured .msh files
  /**
   * Reads .msh files version 4.x and constructs a grid from the cells stored in the file
   * Additionally, stored data can be read.
   *
   * \tparam Grid         Type of the grid to construct.
   * \tparam GridCreator  Policy type to control what to pass to a grid factory with
   *                      data given from the file. See \ref GridCreatorInterface
   *                      [ContinuousGridCreator]
   * \tparam SizeType     Data type for indices in the .msh file. Should match in size
   *                      the `data_size` information in the header of the .msh file.
   *                      This type is used for the internal data in the caches.
   *                      [std::size_t]
   **/
  template <class Grid, class GridCreator = ContinuousGridCreator<Grid>, class SizeType = std::size_t>
  class Gmsh4Reader
    : public FileReader<Grid, Gmsh4Reader<Grid, GridCreator> >
  {
    // Sections visited during the file parsing
    enum class Sections {
      NO_SECTION = 0, MESH_FORMAT, PHYSICAL_NAMES, ENTITIES, PARTITIONED_ENTITIES, NODES, ELEMENTS, PERIODIC,
      GHOST_ELEMENTS, PARAMETRIZATION, NODE_DATA, ELEMENT_DATA, ELEMENT_NODE_DATA, INTERPOLATION_SCHEME
    };

    using Entity = typename Grid::template Codim<0>::Entity;
    using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;

  public:
    using size_type = SizeType;

    /// Constructor. Creates a new GridCreator with the passed factory
    template <class ... Args,
        std::enable_if_t<std::is_constructible<GridCreator, Args...>::value,int> = 0>
    explicit Gmsh4Reader (Args&&... args)
      : creator_(std::make_shared<GridCreator>(std::forward<Args>(args)...))
    {}

    /// Constructor. Stores the references in a non-destroying shared_ptr
    explicit Gmsh4Reader (GridCreator& creator)
      : creator_(stackobject_to_shared_ptr(creator))
    {}

    /// Constructor. Stores the shared_ptr
    explicit Gmsh4Reader (std::shared_ptr<GridCreator> creator)
      : creator_(std::move(creator))
    {}

    /// Read the grid from file with `filename` into the GridFactory stored in the GridCreator
    /**
     * \param filename       The name of the input file
     * \param fillCreator    If `false`, only fill internal data structures,
     *                       if `true`, also calls \ref fillGridCreator. [true]
     **/
    void read (std::string const& filename, bool fillCreator = true);

    /// Advanced read methods
    /// @{

    /// Read the grid from an input stream, referring to a .msh file, into the GridCreator
    /**
     * \param input          A STL input stream to read the Gmsh file from.
     * \param fillCreator    If `false`, only fill internal data structures,
     *                       if `true`, also calls \ref fillGridCreator. [true]
     **/
    void readSerialFileFromStream (std::ifstream& input, bool fillCreator = true);

    /// Read the grid from and input stream, referring to a .pro file, into the GridCreator
    /**
     * \param input          A STL input stream to read the Gmsh file from.
     * \param fillCreator    If `false`, only fill internal data structures,
     *                       if `true`, also calls \ref fillGridCreator. [true]
     **/
    void readParallelFileFromStream (std::ifstream& input, int rank, int size, bool fillCreator = true);

    /// Construct a grid using the GridCreator
    // NOTE: requires the internal data structures to be filled by an aforgoing call to read
    void fillGridCreator (bool insertPieces = true);

    /// @}

    /// Return the filenames of parallel pieces
    std::vector<std::string> const& pieces () const
    {
      return pieces_;
    }

#ifndef DOXYGEN
    /// Implementation of \ref FileReader interface
    static void fillFactoryImpl (GridFactory<Grid>& factory, std::string const& filename)
    {
      Gmsh4Reader reader{factory};
      reader.read(filename);
    }
#endif

  private:
    template<class T>
    void readValueBinary(std::ifstream& input, T &v);
    void readMeshFormat (std::ifstream& input,
                         double& version, int& file_type, int& data_size);

    void readPhysicalNames (std::ifstream& input);
    void readEntitiesAscii (std::ifstream& input);
    void readEntitiesBinary (std::ifstream& input);
    void readPartitionedEntitiesAscii (std::ifstream& input);
    void readPartitionedEntitiesBinary (std::ifstream& input);
    void readNodesAscii (std::ifstream& input);
    void readNodesBinary (std::ifstream& input);
    void readElementsAscii (std::ifstream& input);
    void readElementsBinary (std::ifstream& input);
    void readPeriodic (std::ifstream& input);
    void readGhostElements (std::ifstream& input);
    void readParametrization (std::ifstream& input);
    void readNodeData (std::ifstream& input);
    void readElementData (std::ifstream& input);
    void readElementNodeData (std::ifstream& input);
    void readInterpolationScheme (std::ifstream& input);

    // Test whether line belongs to section
    bool isSection (std::string line,
                    std::string key,
                    Sections current,
                    Sections parent = Sections::NO_SECTION) const
    {
      bool result = line.substr(1, key.length()) == key;
      if (result && current != parent)
        DUNE_THROW(Dune::Exception , "<" << key << "> in wrong section." );
      return result;
    }

    void readString(std::istream& /*stream*/, std::string& name)
    {
      DUNE_THROW(Dune::NotImplemented, "Method readString() not yet implemented");
      name = "";
    }

    // clear all vectors
    void clear ()
    {
      pieces_.clear();

      physicalNames_.clear();

      points_.clear();
      curves_.clear();
      surfaces_.clear();
      volumes_.clear();

      numPartitions_ = 0;
      ghostEntities_.clear();
      partitionedPoints_.clear();
      partitionedCurves_.clear();
      partitionedSurfaces_.clear();
      partitionedVolumes_.clear();

      numNodes_ = 0;
      minNodeTag_ = 0;
      maxNodeTag_ = 0;
      nodes_.clear();

      numElements_ = 0;
      minElementTag_ = 0;
      maxElementTag_ = 0;
      elements_.clear();

      periodic_.clear();
      ghostElements_.clear();
      parametrization_.clear();
      nodeData_.clear();
      elementData_.clear();
      elementNodeData_.clear();
      interpolationScheme_.clear();
    }

    auto comm () const
    {
      return MPIHelper::getCommunication();
    }

  private: // structures used to store data from file

    struct PhysicalNamesAttributes;
    struct PointAttributes;
    struct EntityAttributes;

    template <class Attr>
    struct PartitionedAttributes : public Attr
    {
      int parentDim;
      int parentTag;
      std::vector<int> partitions;
    };

    struct GhostAttributes;
    struct NodeAttributes;
    struct ElementAttributes;
    struct PeriodicAttributes;

    struct GhostElementAttributes {};
    struct ParametrizationAttributes {};
    struct NodeDataAttributes {};
    struct ElementDataAttributes {};
    struct ElementNodeDataAttributes {};
    struct InterpolationSchemeAttributes {};

  private:
    std::shared_ptr<GridCreator> creator_ = nullptr;

    // swap will be true if a binary msh-file has different endianness as the user's system.
    bool swap = false;

    std::vector<std::string> pieces_;

    // PhysicalNames section
    std::vector<PhysicalNamesAttributes> physicalNames_;

    // Entities section
    std::vector<PointAttributes> points_;
    std::vector<EntityAttributes> curves_;
    std::vector<EntityAttributes> surfaces_;
    std::vector<EntityAttributes> volumes_;

    // PartitionedEntities section
    size_type numPartitions_ = 0;
    std::vector<GhostAttributes> ghostEntities_;
    std::vector<PartitionedAttributes<PointAttributes> > partitionedPoints_;
    std::vector<PartitionedAttributes<EntityAttributes> > partitionedCurves_;
    std::vector<PartitionedAttributes<EntityAttributes> > partitionedSurfaces_;
    std::vector<PartitionedAttributes<EntityAttributes> > partitionedVolumes_;

    size_type numNodes_ = 0;
    size_type minNodeTag_ = 0;
    size_type maxNodeTag_ = 0;
    std::vector<NodeAttributes> nodes_;

    size_type numElements_ = 0;
    size_type minElementTag_ = 0;
    size_type maxElementTag_ = 0;
    std::vector<ElementAttributes> elements_;
    std::vector<PeriodicAttributes> periodic_;
    std::vector<GhostElementAttributes> ghostElements_;
    std::vector<ParametrizationAttributes> parametrization_;
    std::vector<NodeDataAttributes> nodeData_;
    std::vector<ElementDataAttributes> elementData_;
    std::vector<ElementNodeDataAttributes> elementNodeData_;
    std::vector<InterpolationSchemeAttributes> interpolationScheme_;

    /// Map the elementType number to number of nodes
    static std::map<int, size_type> elementType_;

    // Associate a string with the corresponding Sections enum
    static std::map<std::string, Sections> sections_;
  };

  // deduction guides
  template <class Grid>
  Gmsh4Reader (GridFactory<Grid>&)
  -> Gmsh4Reader<Grid, ContinuousGridCreator<Grid> >;

  template <class GridCreator,
      class = std::void_t<typename GridCreator::Grid> >
  Gmsh4Reader (GridCreator&)
  -> Gmsh4Reader<typename GridCreator::Grid, GridCreator>;

  template <class GridCreator,
      class = std::void_t<typename GridCreator::Grid> >
  Gmsh4Reader (std::shared_ptr<GridCreator>)
  -> Gmsh4Reader<typename GridCreator::Grid, GridCreator>;

} // end namespace Dune::Impl::Gmsh

#include "gmsh4reader.impl.hh"

#endif
