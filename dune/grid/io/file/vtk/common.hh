// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_COMMON_HH
#define DUNE_GRID_IO_FILE_VTK_COMMON_HH

#include <limits>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
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
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum OutputType {
      //! Output to the file is in ascii.
      ascii,
      //! Output to the file is inline base64 binary.
      base64,
      //! Output is to the file is appended raw binary
      appendedraw,
      //! Output is to the file is appended base64 binary
      appendedbase64
      // //! Output to the file is compressed inline binary.
      // binarycompressed,
      // //! Output is compressed and appended to the file.
      // compressedappended
    };
    //! Whether to produce conforming or non-conforming output.
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
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
       * Each element has its own set of vertices.  The position of a vertex
       * of one element will coincide with the position of the corresponding
       * vertex on another element.  This allows for multiple DoFs (one per
       * element) on the "same" vertex.
       */
      nonconforming
    };

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
      typedef std::conditional<std::numeric_limits<char>::is_signed,
          int, unsigned>::type
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
     * #include <dune/grid/io/file/vtk/common.hh>
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
      polygon = 7,
      quadrilateral = 9,
      tetrahedron = 10,
      hexahedron = 12,
      prism = 13,
      pyramid = 14,
      polyhedron = 42
    };

    //! mapping from GeometryType to VTKGeometryType
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
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

      if (t.isNone() )
      {
        if( t.dim() == 2 ) return polygon;
        if( t.dim() == 3 ) return polyhedron;
      }

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
     * #include <dune/grid/io/file/vtk/common.hh>
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
     * #include <dune/grid/io/file/vtk/common.hh>
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
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum FileType {
      //! for .vtp files (PolyData)
      polyData,
      //! for .vtu files (UnstructuredGrid)
      unstructuredGrid
    };


    //! Descriptor struct for VTK fields
    /**
     * This struct provides general information about a data field to be
     * written to a VTK file.
     *
     * It currently stores the data type and the number of components as well as
     * the name of the field.
     */
    class FieldInfo
    {

    public:

      //! VTK data type
      enum class Type {
        //! scalar field (may also be multi-component, but is treated as a simply
        //! array by ParaView
        scalar,
        //! vector-valued field (always 3D, will be padded if necessary)
        vector,
        //! tensor field (always 3x3)
        tensor
      };

      //! Create a FieldInfo instance with the given name, type and size.
      FieldInfo(std::string name, Type type, std::size_t size)
        : _name(name)
        , _type(type)
        , _size(size)
      {}

      //! The name of the data field
      std::string name() const
      {
        return _name;
      }

      //! The type of the data field
      Type type() const
      {
        return _type;
      }

      //! The number of components in the data field.
      std::size_t size() const
      {
        return _size;
      }

    private:

      std::string _name;
      Type _type;
      std::size_t _size;

    };


  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_COMMON_HH
