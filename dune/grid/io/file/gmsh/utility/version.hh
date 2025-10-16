// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_UTILITY_VERSION_HH
#define DUNE_GRID_IO_FILE_GMSH_UTILITY_VERSION_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include "string.hh"

namespace Dune::Impl::Gmsh
{
  /// Return a version tuple identifying the .msh file format version
  inline std::vector<int> fileVersion(std::string filename)
  {
    std::ifstream file(filename, std::ios_base::in);
    std::string section;
    file >> section;

    if (section != "$MeshFormat")
      DUNE_THROW(Dune::IOError, "Invalid header of msh file.");

    std::string version;
    int file_type = -1;
    int data_size = -1;
    file >> version >> file_type >> data_size;

    if (std::stod(version) <= 0.0)
      DUNE_THROW(Dune::IOError, "Invalid version number in msh file.");

    if (file_type != 0 and file_type != 1)
      DUNE_THROW(Dune::IOError, "Invalid file-type: 0 for ASCII mode, 1 for binary mode.");

    if (data_size < 4 || data_size > 16)
      DUNE_THROW(Dune::IOError, "Invalid data-size range: should be in {4, 16}");

    std::vector<int> version_tuple;
    split(version.begin(), version.end(), '.', [&](auto first, auto last) {
      version_tuple.push_back(std::stoi(std::string{first,last}));
    });

    return version_tuple;
  }
}

#endif
