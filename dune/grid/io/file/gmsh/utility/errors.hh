// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_GRID_IO_FILE_GMSH_UTILITY_ERRORS_HH
#define DUNE_GRID_IO_FILE_GMSH_UTILITY_ERRORS_HH

#include <dune/common/exceptions.hh>

/**
 * \file
 * \brief Macro for wrapping error checks and throwing exceptions
 */

namespace Dune::Impl::Gmsh
{

  class Gmsh4Error : public Exception {};

}

/**
 * \brief check if condition \a cond holds; otherwise, throw a Gmsh4Error with a message.
 */
#define GMSH4_ASSERT_MSG(cond, text)    \
        do {                                  \
          if (!(cond))                        \
          DUNE_THROW(Dune::Impl::Gmsh::Gmsh4Error, text); \
        } while (false)


/**
 * \brief check if condition \a cond holds; otherwise, throw a Gmsh4Error.
 */
#define GMSH4_ASSERT(cond)               \
        do {                                   \
          if (!(cond))                         \
          DUNE_THROW(Dune::Impl::Gmsh::Gmsh4Error, #cond); \
        } while (false)

#endif
