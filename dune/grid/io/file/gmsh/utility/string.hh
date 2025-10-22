// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_UTILITY_STRING_HH
#define DUNE_GRID_IO_FILE_GMSH_UTILITY_STRING_HH

#include <algorithm>
#include <cctype>
#include <locale>
#include <sstream>
#include <string>

namespace Dune::Impl::Gmsh
{
  /// trim a string from the left
  inline std::string& ltrim(std::string& str)
  {
    auto it =  std::find_if(str.begin(), str.end(), [](char ch)
    {
      return !std::isspace<char>(ch, std::locale::classic());
    });
    str.erase(str.begin() , it);
    return str;
  }

  /// trim a string from the right
  inline std::string& rtrim(std::string& str)
  {
    auto it =  std::find_if(str.rbegin(), str.rend(), [](char ch)
    {
      return !std::isspace<char>(ch, std::locale::classic());
    });
    str.erase(it.base(), str.end());
    return str;
  }

  /// trim a string from both sides
  inline std::string& trim(std::string& str)
  {
    return ltrim(rtrim(str));
  }

  template <class InputIter, class T, class Func>
  void split(InputIter first, InputIter end, T const& t, Func f)
  {
    if (first == end)
      return;

    while (true) {
      InputIter found = std::find(first, end, t);
      f(first, found);
      if (found == end)
        break;
      first = ++found;
    }
  }

} // end namespace Dune::Impl::Gmsh

#endif
