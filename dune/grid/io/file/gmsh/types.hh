// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_TYPES_HH
#define DUNE_GRID_IO_FILE_GMSH_TYPES_HH

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <dune/common/ftraits.hh>
#include <dune/geometry/type.hh>

namespace Dune::Impl::Gmsh
{
  /// Get the Dune GeometryType corresponding to a Gmsh entity number
  GeometryType gmshNumberToGeometryType (int elementType);

  /// Mapping of Dune geometry types to Gmsh cell types
  class CellType
  {
  public:
    CellType (GeometryType const& t);

    /// Return Gmsh Cell type
    GeometryType type () const
    {
      return type_;
    }

    /// Return the Dune local number of a cell's vertex given in Gmsh numbering
    int gmshVertexToDuneVertex (int idx) const
    {
      return permutation_[idx];
    }

    /// Return true if Dune and Gmsh use the same vertex numbering for this element type
    bool identityPermutation () const
    {
      return identityPermutation_;
    }

  private:
    GeometryType type_;
    std::vector<int> permutation_;
    bool identityPermutation_;
  };

} // end namespace Dune::Impl::Gmsh

#endif
