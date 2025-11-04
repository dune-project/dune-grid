// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_GMSHREADER_HH
#define DUNE_GRID_IO_FILE_GMSHREADER_HH

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmsh/gmsh2parser.hh>
#include <dune/grid/io/file/gmsh/gmsh4reader.hh>
#include <dune/grid/io/file/gmsh/utility/version.hh>

namespace Dune
{

  /**
     \ingroup Gmsh
     \{
   */

  /** \brief Options for read operation
   *
   * \deprecated The class is deprecated, because it looks like it should be used
   * in the GmshReader implementation, but it is not actually used anywhere.
   */
  struct [[deprecated]] GmshReaderOptions
  {
    enum GeometryOrder {
      /** @brief edges are straight lines. */
      firstOrder,
      /** @brief quadratic boundary approximation. */
      secondOrder
    };
  };

  namespace Gmsh {
    /**
      \ingroup Gmsh
      \brief Options for the Gmsh mesh file reader
      \note Two or more options can be composed using the binary "|"-operator
    **/
    enum class ReaderOptions
    {
      verbose = 1,
      insertBoundarySegments = 2,
      readElementData = 4,
      readBoundaryData = 8
    };

    //! composition operator for reader options
    constexpr ReaderOptions operator | (ReaderOptions a, ReaderOptions b)
    {
      return static_cast<ReaderOptions>(
        static_cast<int>(a) | static_cast<int>(b)
      );
    }

    //! query operator for reader options (is b set in a)
    constexpr bool operator & (ReaderOptions a, ReaderOptions b)
    {
      return static_cast<int>(a) & static_cast<int>(b);
    }

  } // end namespace Gmsh

  /** \brief The GmshReaderParser class has been renamed and moved to the Impl::Gmsh namespace
   *
   * With its new location, the class is now officially not part of the public interface.
   * Please complain if you disagree.
   *
   * \deprecated The name GmshReaderParser remains for backward compatibility,
   * and will be removed eventually.
   *
   * \since dune-grid 2.11
   */
  template <typename GridType>
  struct [[deprecated("Public interface of the GmshReaderParser has been deprecated since dune 2.11.")]]
  GmshReaderParser : Impl::Gmsh::Gmsh2Parser<GridType>
  {
    using Impl::Gmsh::Gmsh2Parser<GridType>::Gmsh2Parser;
  };

  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file

     Read a .msh (version 2) file generated using Gmsh and construct a grid using the grid factory interface.

     The file format used by gmsh can hold grids that are more general than the simplex grids that
     the gmsh grid generator is able to construct.  We try to read as many grids as possible, as
     long as they are valid files.  You can test this by checking whether gmsh will load the file
     and display its content.

     All grids in a gmsh file live in three-dimensional Euclidean space.  If the world dimension
     of the grid type that you are reading the file into is less than three, the remaining coordinates
     are simply ignored.
   */
  template<typename GridType>
  class GmshReader
  {
    //! internal general reading method
    /**
     * This method does all the highlevel steering of the reader:
     * - it will register the GmshReader boundary segment implementation with
     *   the factory
     * - it will ensure the reader is called on all ranks (debug mode only)
     * - proceed to construct a parser (rank 0 only)
     * - use the parser to read the grid into the factory (rank 0 only)
     * - move entity and boundary data from the parser into the data vector
     *   arguments, or clear the data vector arguments, depending on rank
     *
     * \note That the parser always reads the data vectors from the files.
     *       However, if insertBoundarySegments is false, no boundary segments
     *       are inserted into the factory, and thus there will be no correct
     *       indexing of the boundarySegmentToPhysicalEntity vector possible.
     *       For this reason, this method is not exposed to the user, and the
     *       interface methods are responsible to ensure that
     *       boundarySegmentToPhysicalEntity is discarded if boundary segments
     *       are not inserted.
     */
    static void doRead(Dune::GridFactory<GridType> &factory,
                       const std::string &fileName,
                       std::vector<int>& boundarySegmentToPhysicalEntity,
                       std::vector<int>& elementToPhysicalEntity,
                       bool verbose, bool insertBoundarySegments)
    {
      // register boundary segment to boundary segment factory for possible load balancing
      // this needs to be done on all cores since the type might not be known otherwise
      Impl::Gmsh::GmshReaderQuadraticBoundarySegment< Grid::dimension, Grid::dimensionworld >::registerFactory();

#ifndef NDEBUG
      // check that this method is called on all cores
      factory.comm().barrier();
#endif

      // create parse object and read grid on process 0
      if (factory.comm().rank() == 0)
      {
        Impl::Gmsh::Gmsh2Parser<Grid> parser(factory,verbose,insertBoundarySegments);
        parser.read(fileName);

        boundarySegmentToPhysicalEntity = std::move(parser.boundaryIdMap());
        elementToPhysicalEntity = std::move(parser.elementIndexMap());
      }
      else
      {
        boundarySegmentToPhysicalEntity = {};
        elementToPhysicalEntity = {};
      }
    }

    //! pass a discarded lvalue argument to a function
    /**
     * This method is intended to be used in function calls that require
     * lvalue arguments, when the caller just wants to pass in temporary
     * variable that is immediately discarded after the return of the
     * function.  It expects an rvalue argument, that is turned into an
     * lvalue.  For instance:
     * ```c++
     * doRead(factory, fileName, discarded(std::vector<int>{}),
     *        discarded(std::vector<int>{}));
     * ```
     * Here, the vectors are constructed as rvalues, passed through
     * `discarded()` which turns them into lvalues, so they can be arguments
     * to `doRead()`.  `doRead()` will fill them with some data, and they
     * will be destroyed at the end of the full-expression containing the
     * function call.
     *
     * \note It is very likely an error to use this outside a function call
     *       argument.
     */
    template<class T>
    static T &discarded(T &&value) { return static_cast<T&>(value); }

    struct DataArg {
      std::vector<int> *data_ = nullptr;
      DataArg(std::vector<int> &data) : data_(&data) {}
      DataArg(const decltype(std::ignore)&) {}
      DataArg() = default;
    };

    struct DataFlagArg : DataArg {
      bool flag_ = false;
      using DataArg::DataArg;
      DataFlagArg(bool flag) : flag_(flag) {}
    };

  public:
    typedef GridType Grid;

    /** \todo doc me
     */
    static std::unique_ptr<Grid> read (const std::string& fileName, bool verbose = true, bool insertBoundarySegments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      read(factory, fileName, verbose, insertBoundarySegments);

      return factory.createGrid();
    }

    /**
     * \brief Read Gmsh file, possibly with data
     * \param fileName                        Name of the file to read from.
     * \param boundarySegmentToPhysicalEntity Container to fill with boundary segment
     *                                        physical entity data (if insertBoundarySegments=true)
     * \param elementToPhysicalEntity         Container to fill with element physical entity data
     * \param verbose                         Whether to be chatty
     * \param insertBoundarySegments          Whether boundary segments are inserted into the factory
     *
     * \note When insertBoundarySegments=false there is no way to correctly use the values returned
     *       in boundarySegmentToPhysicalEntity. Make sure to set insertBoundarySegments=true if you
     *       intent to do this. An alternative is to use the other overloads which provide compile-time
     *       checking of the provided parameter combinations.
     *
     * \warning The following does not work when the file is in Version 4 Gmsh format:
     * Filling the boundarySegmentToPhysicalEntity and elementToPhysicalEntity fields,
     * controlling verbosity, and inserting boundary segments.
     *
     * \todo This interface is error-prone and should not be exposed to the user. However, the
     *       compile-time overloads may not provide sufficient runtime flexibility in all cases.
     *       Therefore this interface is kept until a better interface can be agreed on.
     *       See https://gitlab.dune-project.org/core/dune-grid/-/issues/107
     */
    static std::unique_ptr<Grid> read (const std::string& fileName,
                       std::vector<int>& boundarySegmentToPhysicalEntity,
                       std::vector<int>& elementToPhysicalEntity,
                       bool verbose = true, bool insertBoundarySegments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      if (Impl::Gmsh::fileVersion(fileName)[0]==4)
      {
        Impl::Gmsh::Gmsh4Reader<Grid>::fillFactory(factory, fileName);
        return factory.createGrid();
      }

      doRead(
        factory, fileName, boundarySegmentToPhysicalEntity,
        elementToPhysicalEntity, verbose, insertBoundarySegments
      );

      return factory.createGrid();
    }

    /** \brief Read Gmsh grid file into a `GridFactory` object
     *
     * \warning Controlling verbosity, and inserting boundary segments does not work
     * when the file is in Version 4 Gmsh format.
     */
    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
                      bool verbose = true, bool insertBoundarySegments=true)
    {
      if (Impl::Gmsh::fileVersion(fileName)[0]==4)
      {
        Impl::Gmsh::Gmsh4Reader<Grid>::fillFactory(factory, fileName);
        return;
      }

      doRead(
        factory, fileName, discarded(std::vector<int>{}),
        discarded(std::vector<int>{}), verbose, insertBoundarySegments
      );
    }

    //! read Gmsh file, possibly with data
    /**
     * \param factory             The GridFactory to fill.
     * \param fileName            Name of the file to read from.
     * \param boundarySegmentData Container to fill with boundary segment
     *                            physical entity data, or `std::ignore`, or a
     *                            `bool` value.  Boundary segments are
     *                            inserted when a container or `true` is
     *                            given, otherwise they are not inserted.
     * \param elementData         Container to fill with element physical
     *                            entity data, or `std::ignore`.
     * \param verbose             Whether to be chatty.
     *
     * Containers to fill with data must be `std::vector<int>` lvalues.
     * Element data is indexed by the insertion index of the element,
     * boundarySegment data is indexed by the insertion index of the boundary
     * intersection.  These can be obtained from the `factory`, and are lost
     * once the grid gets modified (refined or load-balanced).
     *
     * \warning The following does not work when the file is in Version 4 Gmsh format:
     * Filling the boundarySegmentData and elementData fields, and controlling verbosity.
     *
     * \note At the moment the data containers are still filled internally,
     *       even if they are ignored.  So not having to pass them is more of
     *       a convenience feature and less of an optimization.  This may
     *       however change in the future.
     */
    static void read (Dune::GridFactory<Grid> &factory,
                      const std::string &fileName,
                      DataFlagArg boundarySegmentData,
                      DataArg elementData,
                      bool verbose=true)
    {
      if (Impl::Gmsh::fileVersion(fileName)[0]==4)
      {
        Impl::Gmsh::Gmsh4Reader<Grid>::fillFactory(factory, fileName);
        return;
      }

      doRead(
        factory, fileName,
        boundarySegmentData.data_
          ? *boundarySegmentData.data_ : discarded(std::vector<int>{}),
        elementData.data_
          ? *elementData.data_ : discarded(std::vector<int>{}),
        verbose,
        boundarySegmentData.flag_ || boundarySegmentData.data_
      );
    }

    /**
     * \brief Read Gmsh file, possibly with data
     * \param factory                         The GridFactory to fill.
     * \param fileName                        Name of the file to read from.
     * \param boundarySegmentToPhysicalEntity Container to fill with boundary segment
     *                                        physical entity data (if insertBoundarySegments=true)
     * \param elementToPhysicalEntity         Container to fill with element physical entity data
     * \param verbose                         Whether to be chatty
     * \param insertBoundarySegments          Whether boundary segments are inserted into the factory
     *
     * \note When insertBoundarySegments=false there is no way to correctly use the values returned
     *       in boundarySegmentToPhysicalEntity. Make sure to set insertBoundarySegments=true if you
     *       intent to do this. An alternative is to use the other overloads which provide compile-time
     *       checking of the provided parameter combinations.
     *
     * \warning The following does not work when the file is in Version 4 Gmsh format:
     * Filling the boundarySegmentToPhysicalEntity and elementToPhysicalEntity fields,
     * controlling verbosity, and inserting boundary segments.
     *
     * \todo This interface is error-prone and should not be exposed to the user. However, the
     *       compile-time overloads may not provide sufficient runtime flexibility in all cases.
     *       Therefore this interface is kept until a better interface can be agreed on.
     *       See https://gitlab.dune-project.org/core/dune-grid/-/issues/107
     */
    static void read (Dune::GridFactory<Grid>& factory,
                      const std::string& fileName,
                      std::vector<int>& boundarySegmentToPhysicalEntity,
                      std::vector<int>& elementToPhysicalEntity,
                      bool verbose, bool insertBoundarySegments)
    {
      if (Impl::Gmsh::fileVersion(fileName)[0]==4)
      {
        Impl::Gmsh::Gmsh4Reader<Grid>::fillFactory(factory, fileName);
        return;
      }

      doRead(
        factory, fileName, boundarySegmentToPhysicalEntity,
        elementToPhysicalEntity, verbose, insertBoundarySegments
      );
    }

    //! Dynamic Gmsh reader interface
    //\{

    using Opts = Gmsh::ReaderOptions;

    static constexpr Opts defaultOpts =
      Opts::verbose | Opts::insertBoundarySegments | Opts::readElementData | Opts::readBoundaryData;

    //! Construct a Gmsh reader object (alternatively use one of the static member functions)

    /**
     * \brief Construct a Gmsh reader object from a file name
     * \param fileName Name of the file to read from.
     * \param options Options of the type `Dune::Gmsh::ReaderOptions`
     *
     * To pass several options, combine them with the |-operator like this
     *
      \code
      using Opt = Dune::Gmsh::ReaderOptions;
      auto reader = Dune::GmshReader("grid.msh", Opt::verbose | Opt::readElementData)
      \endcode
     *
     * Per default the reader has enabled the following options
       - Dune::Gmsh::ReaderOptions::verbose
       - Dune::Gmsh::ReaderOptions::insertBoundarySegments
       - Dune::Gmsh::ReaderOptions::readBoundaryData
       - Dune::Gmsh::ReaderOptions::readElementData
     *
     * Passing any option to the interface will overwrite these defaults.
     *
     * A Dune grid object can be obtained via the `createGrid()` member
     *
     * \warning All options are automatically `false` when the file is in Version 4 Gmsh format.
     */
    GmshReader(const std::string& fileName,
               Gmsh::ReaderOptions options = defaultOpts)
    {
      gridFactory_ = std::make_unique<Dune::GridFactory<Grid>>();
      readGridFile(fileName, *gridFactory_, options);
    }

    /**
     * \brief Construct a Gmsh reader object from a file name and a grid factory
     * \param fileName Name of the file to read from.
     * \param options Options of the type `Dune::Gmsh::ReaderOptions`
     *
     * Use this constructor if you need access to the grid factory after the grid
     * has been read, e.g., for obtaining boundary segment insertion indices.
     */
    GmshReader(const std::string& fileName, GridFactory<Grid>& factory,
               Gmsh::ReaderOptions options = defaultOpts)
    {
      readGridFile(fileName, factory, options);
    }

    //! Access element data (maps element index to Gmsh physical entity)
    const std::vector<int>& elementData () const
    {
      checkElementData();
      return elementIndexToGmshPhysicalEntity_;
    }

    //! Access boundary data (maps boundary segment index to Gmsh physical entity)
    const std::vector<int>& boundaryData () const
    {
      checkBoundaryData();
      return boundarySegmentIndexToGmshPhysicalEntity_;
    }

    /**
     * \brief If element data is available
     * \note This is false if no such data was requested
     */
    bool hasElementData () const
    { return hasElementData_ && !extractedElementData_; }

    /**
     * \brief If boundary data is available
     * \note This is false if no such data was requested
     */
    bool hasBoundaryData () const
    { return hasBoundaryData_ && !extractedBoundaryData_; }

    //! Erase element data from reader and return the data
    std::vector<int> extractElementData ()
    {
      checkElementData();
      extractedElementData_ = true;
      return std::move(elementIndexToGmshPhysicalEntity_);
    }

    //! Erase boundary data from reader and return the data
    std::vector<int> extractBoundaryData ()
    {
      checkBoundaryData();
      extractedBoundaryData_ = true;
      return std::move(boundarySegmentIndexToGmshPhysicalEntity_);
    }

    //! Create the grid
    std::unique_ptr<Grid> createGrid ()
    {
      if (!gridFactory_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed with a Dune::GridFactory. "
          << "This grid factory has been filled with all information to create a grid. "
          << "Please use this factory to create the grid by calling factory.createGrid(). "
          << "Alternatively use the constructor without passing the factory in combination with this member function."
        );

      return gridFactory_->createGrid();
    }

    //\}

  private:
    void checkElementData () const
    {
      if (!hasElementData_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed without the option 'readElementData'. "
          << "Please enable reading element data by passing the option 'Gmsh::ReaderOpts::readElementData' "
          << "to the constructor of this class."
        );

      if (extractedElementData_)
        DUNE_THROW(Dune::InvalidStateException,
          "The element data has already been extracted from this GmshReader "
          << "via a function call to reader.extractElementData(). Use the extracted data or "
          << "read the grid data from file again by constructing a new reader."
        );
    }

    void checkBoundaryData () const
    {
      if (!hasBoundaryData_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed without the option 'readBoundaryData'. "
          << "Please enable reading boundary data by passing the option 'Gmsh::ReaderOpts::readBoundaryData' "
          << "to the constructor of this class."
        );

      if (extractedBoundaryData_)
        DUNE_THROW(Dune::InvalidStateException,
          "The boundary data has already been extracted from this GmshReader "
          << "via a function call to reader.extractBoundaryData(). Use the extracted data or "
          << "read the grid data from file again by constructing a new reader."
        );
    }

    void readGridFile (const std::string& fileName, GridFactory<Grid>& factory, Gmsh::ReaderOptions options)
    {
      if (Impl::Gmsh::fileVersion(fileName)[0]==4)
      {
        Impl::Gmsh::Gmsh4Reader<Grid>::fillFactory(factory, fileName);
        return;
      }

      const bool verbose = options & Opts::verbose;
      const bool insertBoundarySegments = options & Opts::insertBoundarySegments;
      const bool readBoundaryData = options & Opts::readBoundaryData;
      const bool readElementData = options & Opts::readElementData;

      doRead(
        factory, fileName, boundarySegmentIndexToGmshPhysicalEntity_,
        elementIndexToGmshPhysicalEntity_, verbose,
        readBoundaryData || insertBoundarySegments
      );

      // clear unwanted data
      if (!readBoundaryData)
          boundarySegmentIndexToGmshPhysicalEntity_ = std::vector<int>{};
      if (!readElementData)
          elementIndexToGmshPhysicalEntity_ = std::vector<int>{};

      hasElementData_ = readElementData;
      hasBoundaryData_ = readBoundaryData;
    }

    std::unique_ptr<Dune::GridFactory<Grid>> gridFactory_;

    std::vector<int> elementIndexToGmshPhysicalEntity_;
    std::vector<int> boundarySegmentIndexToGmshPhysicalEntity_;

    bool hasElementData_ = false;
    bool hasBoundaryData_ = false;

    // for better error messages, we keep track of these separately
    bool extractedElementData_ = false;
    bool extractedBoundaryData_ = false;
  };

  /** \} */

} // namespace Dune

#endif
