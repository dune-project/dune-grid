// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_COMMON_HH
#define DUNE_GRID_IO_FILE_VTK_COMMON_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Common stuff for the VTKWriter

    This file contains common stuff for all instances of VTKWriter.
 */

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  //////////////////////////////////////////////////////////////////////
  //
  //  VTKOptions
  //

  //! options for VTK output
  struct VTKOptions
  {
    //! How the bulk data should be stored in the file
    enum OutputType {
      /** @brief Output to the file is in ascii. */
      ascii,
      /** @brief Output to the file is inline binary. */
      binary,
      /** @brief Ouput is appended binary to the file. */
      binaryappended
      // /** @brief Output to the file is compressed inline binary. */
      // binarycompressed,
      // /** @brief Ouput is compressed and appended to the file. */
      // compressedappended
    };
    //! Whether to produce conforming or non-conforming output.
    /**
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
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  VTKTypeNameTraits
  //

  //! map type to its VTK name in data array
  /**
   * \tparam T The type whose VTK name is requested
   */
  template<class T>
  struct VTKTypeNameTraits {
    //! return VTK name of the type
    /**
     * If the type is not known to VTK, return empty string.
     */
    std::string operator () (){
      return "";
    }
  };

  template<>
  struct VTKTypeNameTraits<char> {
    std::string operator () () {
      return "Int8";
    }
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned char> {
    std::string operator () () {
      return "UInt8";
    }
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<short> {
    std::string operator () () {
      return "Int16";
    }
    typedef short PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned short> {
    std::string operator () () {
      return "UInt16";
    }
    typedef unsigned short PrintType;
  };

  template<>
  struct VTKTypeNameTraits<int> {
    std::string operator () () {
      return "Int32";
    }
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned int> {
    std::string operator () () {
      return "UInt32";
    }
    typedef unsigned int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<float> {
    std::string operator () () {
      return "Float32";
    }
    typedef float PrintType;
  };

  template<>
  struct VTKTypeNameTraits<double> {
    std::string operator () () {
      return "Float64";
    }
    typedef double PrintType;
  };

  //////////////////////////////////////////////////////////////////////
  //
  //  VTKGeometryType related stuff
  //

  //! Type representing VTK's entity geometry types
  /**
   * Only the types which have a corresponding Dune::GeometryType have been
   * included here.  Dune-type names have been used, this mainly makes a
   * difference for vtkPrism, which is known by VTK as VTK_WEDGE.
   */
  enum VTKGeometryType {
    vtkVertex = 1,
    vtkLine = 3,
    vtkTriangle = 5,
    vtkQuadrilateral = 9,
    vtkTetrahedron = 10,
    vtkHexahedron = 12,
    vtkPrism = 13,
    vtkPyramid = 14
  };

  //! mapping from GeometryType to VTKGeometryType
  inline VTKGeometryType vtkType(const GeometryType & t)
  {
    if (t.isVertex()) return vtkVertex;
    if (t.isLine()) return vtkLine;
    if (t.isTriangle()) return vtkTriangle;
    if (t.isQuadrilateral()) return vtkQuadrilateral;
    if (t.isTetrahedron()) return vtkTetrahedron;
    if (t.isPyramid()) return vtkPyramid;
    if (t.isPrism()) return vtkPrism;
    if (t.isHexahedron()) return vtkHexahedron;

    DUNE_THROW(IOError,"VTKWriter: unsupported GeometryType " << t);
  }

  ////////////////////////////////////////////////////////////////////////
  //
  //  Functions for transforming the index of a corner inside an entity
  //  between Dune and VTK
  //

  //! renumber VTK <-> Dune
  /**
   * Since the renumbering never does anything more complex than exchanging
   * two indices, this method works both ways.
   */
  inline int vtkRenumber(const GeometryType &t, int i)
  {
    static const int quadRenumbering[4] = {0,1,3,2};
    static const int cubeRenumbering[8] = {0,1,3,2,4,5,7,6};
    static const int prismRenumbering[6] = {0,2,1,3,5,4};
    static const int pyramidRenumbering[5] = {0,1,3,2,4};
    if (t.isQuadrilateral())
      return quadRenumbering[i];
    if (t.isPyramid())
      return pyramidRenumbering[i];
    if (t.isPrism())
      return prismRenumbering[i];
    if (t.isHexahedron())
      return cubeRenumbering[i];
    return i;
  }
  //! renumber VTK <-> Dune
  /**
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
  int vtkRenumber(const T& t, int i)
  {
    return vtkRenumber(t.type(), i);
  }

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_COMMON_HH
