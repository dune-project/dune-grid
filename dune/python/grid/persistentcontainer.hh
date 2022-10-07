// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_PERSISTENTCONTAINER_HH
#define DUNE_PYTHON_GRID_PERSISTENTCONTAINER_HH

#include <array>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>


#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/python/pybind11/pybind11.h>

#include <iostream>

namespace Dune
{

  namespace Python
  {

    // registerPersistentContainer
    // ---------------------------

    template< class PersistentContainer, class... options >
    inline static void registerPersistentContainer ( pybind11::handle scope, pybind11::class_< PersistentContainer, options... > cls )
    {
      typedef typename PersistentContainer::Grid Grid;
      cls.def( pybind11::init( [] ( Grid &grid, int codim ) { return new PersistentContainer( grid, codim ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "size", [] ( const PersistentContainer &self ) -> int {
          return self.size();
        },
        R"doc(
          return the size of the given persistent container
        )doc" );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_PERSISTENTCONTAINER_HH
