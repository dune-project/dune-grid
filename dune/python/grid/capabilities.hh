// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_CAPABILITIES_HH
#define DUNE_PYTHON_GRID_CAPABILITIES_HH

#include <type_traits>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, class Coordinates >
  class YaspGrid;



  namespace Python
  {

    namespace Capabilities
    {

      using namespace Dune::Capabilities;



      // HasGridFactory
      // --------------

      template< class Grid >
      struct HasGridFactory
        : public std::integral_constant< bool, !isCartesian< Grid >::v >
      {};



      // HasStructuredGridFactory
      // ------------------------

      template< class Grid >
      struct HasStructuredGridFactory
        : public std::integral_constant< bool, HasGridFactory< Grid >::value && (Grid::dimension == Grid::dimensionworld) >
      {};

      template< int dim, class Coordinates >
      struct HasStructuredGridFactory< YaspGrid< dim, Coordinates > >
        : public std::true_type
      {};



      // canIterate
      // ----------

      template< class Grid, int codim >
      struct canIterate
        : public std::integral_constant< bool, hasEntityIterator< Grid, codim >::v >
      {};

    } // namespace Capabilities

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_CAPABILITIES_HH
