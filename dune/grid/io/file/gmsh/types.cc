// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <dune/grid/io/file/gmsh/types.hh>

#include <iostream>

namespace Dune::Impl::Gmsh
{
  GeometryType gmshNumberToGeometryType (int elementType)
  {
    switch (elementType) {
    case 1: return GeometryTypes::line;
    case 2: return GeometryTypes::triangle;
    case 3: return GeometryTypes::quadrilateral;
    case 4: return GeometryTypes::tetrahedron;
    case 5: return GeometryTypes::hexahedron;
    case 6: return GeometryTypes::prism;
    case 7: return GeometryTypes::pyramid;
    case 8: return GeometryTypes::line;
    case 9: return GeometryTypes::triangle;
    case 10: return GeometryTypes::quadrilateral;
    case 11: return GeometryTypes::tetrahedron;
    case 12: return GeometryTypes::hexahedron;
    case 13: return GeometryTypes::prism;
    case 14: return GeometryTypes::pyramid;
    case 15: return GeometryTypes::vertex;
    case 16: return GeometryTypes::quadrilateral;
    case 17: return GeometryTypes::hexahedron;
    case 18: return GeometryTypes::prism;
    case 19: return GeometryTypes::pyramid;
    case 20: return GeometryTypes::triangle;
    case 21: return GeometryTypes::triangle;
    case 22: return GeometryTypes::triangle;
    case 23: return GeometryTypes::triangle;
    case 24: return GeometryTypes::triangle;
    case 25: return GeometryTypes::triangle;
    case 26: return GeometryTypes::line;
    case 27: return GeometryTypes::line;
    case 28: return GeometryTypes::line;
    case 29: return GeometryTypes::tetrahedron;
    case 30: return GeometryTypes::tetrahedron;
    case 31: return GeometryTypes::tetrahedron;
    case 32: return GeometryTypes::tetrahedron;
    case 33: return GeometryTypes::tetrahedron;
    //case 34: return polygon;
    //case 35: return polygon;
    case 36: return GeometryTypes::quadrilateral;
    case 37: return GeometryTypes::quadrilateral;
    case 38: return GeometryTypes::quadrilateral;
    case 39: return GeometryTypes::quadrilateral;
    case 40: return GeometryTypes::quadrilateral;
    case 41: return GeometryTypes::quadrilateral;
    case 42: return GeometryTypes::triangle;
    case 43: return GeometryTypes::triangle;
    case 44: return GeometryTypes::triangle;
    case 45: return GeometryTypes::triangle;
    case 46: return GeometryTypes::triangle;
    case 47: return GeometryTypes::quadrilateral;
    case 48: return GeometryTypes::quadrilateral;
    case 49: return GeometryTypes::quadrilateral;
    case 50: return GeometryTypes::quadrilateral;
    case 51: return GeometryTypes::quadrilateral;
    case 52: return GeometryTypes::triangle;
    case 53: return GeometryTypes::triangle;
    case 54: return GeometryTypes::triangle;
    case 55: return GeometryTypes::triangle;
    case 56: return GeometryTypes::triangle;
    case 57: return GeometryTypes::quadrilateral;
    case 58: return GeometryTypes::quadrilateral;
    case 59: return GeometryTypes::quadrilateral;
    case 60: return GeometryTypes::quadrilateral;
    case 61: return GeometryTypes::quadrilateral;
    case 62: return GeometryTypes::line;
    case 63: return GeometryTypes::line;
    case 64: return GeometryTypes::line;
    case 65: return GeometryTypes::line;
    case 66: return GeometryTypes::line;
    //case 67: return GeometryTypes::line;
    //case 68: return GeometryTypes::triangle;
    //case 69: return polygon;
    //case 70: return line;
    case 71: return GeometryTypes::tetrahedron;
    case 72: return GeometryTypes::tetrahedron;
    case 73: return GeometryTypes::tetrahedron;
    case 74: return GeometryTypes::tetrahedron;
    case 75: return GeometryTypes::tetrahedron;
    case 79: return GeometryTypes::tetrahedron;
    case 80: return GeometryTypes::tetrahedron;
    case 81: return GeometryTypes::tetrahedron;
    case 82: return GeometryTypes::tetrahedron;
    case 83: return GeometryTypes::tetrahedron;
    case 84: return GeometryTypes::line;
    case 85: return GeometryTypes::triangle;
    case 86: return GeometryTypes::quadrilateral;
    case 87: return GeometryTypes::tetrahedron;
    case 88: return GeometryTypes::hexahedron;
    case 89: return GeometryTypes::prism;
    case 90: return GeometryTypes::prism;
    case 91: return GeometryTypes::prism;
    case 92: return GeometryTypes::hexahedron;
    case 93: return GeometryTypes::hexahedron;
    case 94: return GeometryTypes::hexahedron;
    case 95: return GeometryTypes::hexahedron;
    case 96: return GeometryTypes::hexahedron;
    case 97: return GeometryTypes::hexahedron;
    case 98: return GeometryTypes::hexahedron;
    case 99: return GeometryTypes::hexahedron;
    case 100: return GeometryTypes::hexahedron;
    case 101: return GeometryTypes::hexahedron;
    case 102: return GeometryTypes::hexahedron;
    case 103: return GeometryTypes::hexahedron;
    case 104: return GeometryTypes::hexahedron;
    case 105: return GeometryTypes::hexahedron;
    case 106: return GeometryTypes::prism;
    case 107: return GeometryTypes::prism;
    case 108: return GeometryTypes::prism;
    case 109: return GeometryTypes::prism;
    case 110: return GeometryTypes::prism;
    case 111: return GeometryTypes::prism;
    case 112: return GeometryTypes::prism;
    case 113: return GeometryTypes::prism;
    case 114: return GeometryTypes::prism;
    case 115: return GeometryTypes::prism;
    case 116: return GeometryTypes::prism;
    case 117: return GeometryTypes::prism;
    case 118: return GeometryTypes::pyramid;
    case 119: return GeometryTypes::pyramid;
    case 120: return GeometryTypes::pyramid;
    case 121: return GeometryTypes::pyramid;
    case 122: return GeometryTypes::pyramid;
    case 123: return GeometryTypes::pyramid;
    case 124: return GeometryTypes::pyramid;
    case 125: return GeometryTypes::pyramid;
    case 126: return GeometryTypes::pyramid;
    case 127: return GeometryTypes::pyramid;
    case 128: return GeometryTypes::pyramid;
    case 129: return GeometryTypes::pyramid;
    case 130: return GeometryTypes::pyramid;
    case 131: return GeometryTypes::pyramid;
    case 132: return GeometryTypes::pyramid;
    //case 133: return GeometryTypes::vertex;
    //case 134: return GeometryTypes::line;
    //case 135: return GeometryTypes::triangle;
    //case 136: return GeometryTypes::tetrahedron;
    case 137: return GeometryTypes::tetrahedron;
    //case 138: return GeometryTypes::triangle;
    //case 139: return GeometryTypes::tetrahedron;
    case 140: return GeometryTypes::triangle;
    default:
      DUNE_THROW(RangeError, "CellType does not map to GeometryType.");
      std::abort();
    }
  }

  CellType::CellType (GeometryType const& t)
    : identityPermutation_(true)
  {
    if (t.isVertex()) {
      type_ = GeometryTypes::vertex;
      permutation_ = {0};
    }
    else if (t.isLine()) {
      type_ = GeometryTypes::line;
      permutation_ = {0,1};
    }
    else if (t.isTriangle()) {
      type_ = GeometryTypes::triangle;
      permutation_ = {0,1,2};
    }
    else if (t.isQuadrilateral()) {
      type_ = GeometryTypes::quadrilateral;
      permutation_ = {0,1,3,2};
      identityPermutation_ = false;
    }
    else if (t.isTetrahedron()) {
      type_ = GeometryTypes::tetrahedron;
      permutation_ = {0,1,2,3};
    }
    else if (t.isHexahedron()) {
      type_ = GeometryTypes::hexahedron;
      permutation_ = {0,1,3,2,4,5,7,6};
      identityPermutation_ = false;
    }
    else if (t.isPrism()) {
      type_ = GeometryTypes::prism;
      permutation_ = {0,2,1,3,5,4};
      identityPermutation_ = false;
    }
    else if (t.isPyramid()) {
      type_ = GeometryTypes::pyramid;
      permutation_ = {0,1,3,2,4};
      identityPermutation_ = false;
    }
    else if (t.isNone() && t.dim() == 1) {
      type_ = GeometryTypes::line;
      permutation_ = {0,1};
    }
    else {
      std::cerr << "Geometry Type not supported by Gmsh4!\n";
      std::abort();
    }
  }

} // end namespace Dune::Impl::Gmsh
