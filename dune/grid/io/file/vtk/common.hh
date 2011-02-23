// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_COMMON_HH
#define DUNE_GRID_IO_FILE_VTK_COMMON_HH

#include <limits>
#include <sstream>
#include <string>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/typetraits.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Common stuff for the VTKWriter

    This file contains common stuff for all instances of VTKWriter.
 */

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //////////////////////////////////////////////////////////////////////
    //
    //  VTKOptions
    //

    //! How the bulk data should be stored in the file
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum OutputType {
      //! Output to the file is in ascii.
      ascii,
      //! Output to the file is inline base64 binary.
      base64,
      //! Ouput is to the file is appended raw binary
      appendedraw,
      //! Ouput is to the file is appended base64 binary
      appendedbase64
      // //! Output to the file is compressed inline binary.
      // binarycompressed,
      // //! Ouput is compressed and appended to the file.
      // compressedappended
    };
    //! Whether to produce conforming or non-conforming output.
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * This applies to the conformity of the data; a non-conforming grid can
     * still be written in conforming data mode, and it is quite possible for
     * data to be non-conforming on a conforming grid.
     */
    enum DataMode {
      //! Output conforming data.
      /**
       * Neighboring elements share common vertices and thus have a common DoF
       * on that vertex.
       */
      conforming,
      //! Output non-conforming data.
      /**
       * Each element has it's own set of vertices.  The position of a vertex
       * of one element will concide with the position of the corresponding
       * vertex on another element.  This allows for multiple DoFs (one per
       * element) on the "same" vertex.
       */
      nonconforming
    };

  } // namespace VTK

  //! options for VTK output
  /**
   * \code
   *#include <dune/grid/io/file/vtk/common.hh>
   * \endcode
   *
   * \deprecated Use the enums from the namespace VTK instead
   */
  struct DUNE_DEPRECATED VTKOptions {
    //! How the bulk data should be stored in the file
    /**
     * \deprecated Use the enum from the namespace VTK instead.  When
     *             converting your code, look for any occurances of
     *             VTKOptions::ascii, VTKOptions::binary, and
     *             VTKOptions::binaryappended, and convert them es well.
     *             Unfortunately, marking the compatibility constants here as
     *             deprecated seems to have no effect.
     */
    typedef DUNE_DEPRECATED VTK::OutputType OutputType;
    //! Output to the file is in ascii.
    /**
     * \deprecated Use VTK::ascii instead
     */
    static const VTK::OutputType ascii = VTK::ascii;
    //! Output to the file is inline binary.
    /**
     * \deprecated Use VTK::base64 instead
     */
    static const VTK::OutputType binary = VTK::base64;
    //! Ouput is appended binary to the file.
    /**
     * \deprecated Use VTK::appendedraw instead
     */
    static const VTK::OutputType binaryappended = VTK::appendedraw;
    //! Whether to produce conforming or non-conforming output.
    /**
     * This applies to the conformity of the data; a non-conforming grid can
     * still be written in conforming data mode, and it is quite possible for
     * data to be non-conforming on a conforming grid.
     *
     * \deprecated Use the enum from the namespace VTK instead.  When
     *             converting your code, look for any occurances of
     *             VTKOptions::conforming or VTKOptions::nonconforming, and
     *             convert them es well.  Unfortunately, marking the
     *             compatibility constants here as deprecated seems to have no
     *             effect.
     */
    typedef DUNE_DEPRECATED VTK::DataMode DataMode;
    //! Output conforming data.
    /**
     * Neighboring elements share common vertices and thus have a common DoF
     * on that vertex.
     *
     * \deprecated Use the value from the namespace VTK instead
     */
    static const VTK::DataMode conforming = VTK::conforming;
    //! Output non-conforming data.
    /**
     * Each element has it's own set of vertices.  The position of a vertex of
     * one element will concide with the position of the corresponding vertex
     * on another element.  This allows for multiple DoFs (one per element) on
     * the "same" vertex.
     *
     * \deprecated Use the value from the namespace VTK instead
     */
    static const VTK::DataMode nonconforming = VTK::nonconforming;
  };

  namespace VTK {

    //////////////////////////////////////////////////////////////////////
    //
    //  PrintType
    //

    //! determine a type to safely put another type into a stream
    /**
     * This is mainly interating for character types which should print as
     * their integral value, not as a character.
     */
    template<typename T>
    struct PrintType {
      //! type to convert T to before putting it into a stream with <<
      typedef T Type;
    };

    template<>
    struct PrintType<unsigned char> {
      typedef unsigned Type;
    };

    template<>
    struct PrintType<signed char> {
      typedef int Type;
    };

    template<>
    struct PrintType<char> {
      typedef SelectType<std::numeric_limits<char>::is_signed,
          int, unsigned>::Type
      Type;
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  TypeName
    //

    //! map type to its VTK name in data array
    /**
     * \tparam T The type whose VTK name is requested
     */
    template<typename T>
    class TypeName {
      static std::string getString() {
        static const unsigned int_sizes[] = { 8, 16, 32, 64, 0 };
        static const unsigned float_sizes[] = { 32, 64, 0 };
        const unsigned* sizes;

        std::ostringstream s;
        if(std::numeric_limits<T>::is_integer) {
          if(std::numeric_limits<T>::is_signed)
            s << "Int";
          else
            s << "UInt";
          sizes = int_sizes;
        }
        else {
          // assume float
          s << "Float";
          sizes = float_sizes;
        }

        static const unsigned size = 8*sizeof(T);
        while(*sizes != 0 && *sizes <= size) ++sizes;
        --sizes;
        s << *sizes;

        return s.str();
      }

    public:
      //! return VTK name of the type
      /**
       * If the type is not known to VTK, return empty string.
       */
      const std::string& operator()() const {
        static const std::string s = getString();
        return s;
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  VTK::GeometryType related stuff
    //

    //! Type representing VTK's entity geometry types
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * Only the types which have a corresponding Dune::GeometryType have been
     * included here.  Dune-type names have been used, this mainly makes a
     * difference for vtkPrism, which is known by VTK as VTK_WEDGE.
     */
    enum GeometryType {
      vertex = 1,
      line = 3,
      triangle = 5,
      quadrilateral = 9,
      tetrahedron = 10,
      hexahedron = 12,
      prism = 13,
      pyramid = 14
    };

    //! mapping from GeometryType to VTKGeometryType
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    inline GeometryType geometryType(const Dune::GeometryType& t)
    {
      if (t.isVertex()) return vertex;
      if (t.isLine()) return line;
      if (t.isTriangle()) return triangle;
      if (t.isQuadrilateral()) return quadrilateral;
      if (t.isTetrahedron()) return tetrahedron;
      if (t.isPyramid()) return pyramid;
      if (t.isPrism()) return prism;
      if (t.isHexahedron()) return hexahedron;

      DUNE_THROW(IOError,"VTKWriter: unsupported GeometryType " << t);
    }

    ////////////////////////////////////////////////////////////////////////
    //
    //  Functions for transforming the index of a corner inside an entity
    //  between Dune and VTK
    //

    //! renumber VTK <-> Dune
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * Since the renumbering never does anything more complex than exchanging
     * two indices, this method works both ways.
     */
    inline int renumber(const Dune::GeometryType &t, int i)
    {
      static const int quadRenumbering[4] = {0,1,3,2};
      static const int cubeRenumbering[8] = {0,1,3,2,4,5,7,6};
      static const int prismRenumbering[6] = {0,2,1,3,5,4};
      static const int pyramidRenumbering[5] = {0,1,3,2,4};

      if (t.isQuadrilateral()) return quadRenumbering[i];
      if (t.isPyramid()) return pyramidRenumbering[i];
      if (t.isPrism()) return prismRenumbering[i];
      if (t.isHexahedron()) return cubeRenumbering[i];

      return i;
    }

    //! renumber VTK <-> Dune
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * This function is just a convenience shortcut function wrapping
     * renumber(const GeometryType&, int).
     *
     * \param t Entity, Intersection or Geometry to do the renumbering in.
     *          Basically, anything with a method type() returning a
     *          GeometryType should work here.
     * \param i Index to of corner in either Dune or VTK numbering (the result
     *          will be in the other numbering)
     */
    template<typename T>
    int renumber(const T& t, int i)
    {
      return renumber(t.type(), i);
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  Determine Endianness
    //

    //! determine endianness of this C++ implementation
    /**
     * \returns A string suitable for the byte_order property in VTK files;
     *          either "BigEndian" or "LittleEndian".
     */
    inline std::string getEndiannessString()
    {
      short i = 1;
      if (reinterpret_cast<char*>(&i)[1] == 1)
        return "BigEndian";
      else
        return "LittleEndian";
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  which type of vtkfile to write
    //

    //! which type of VTK file to write
    /**
     * \code
     *#include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum FileType {
      //! for .vtp files (PolyData)
      polyData,
      //! for .vtu files (UnstructuredGrid)
      unstructuredGrid
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_COMMON_HH
