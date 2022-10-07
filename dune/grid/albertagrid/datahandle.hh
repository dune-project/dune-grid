// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDDATAHANDLE_HH
#define DUNE_ALBERTAGRIDDATAHANDLE_HH

#include <iostream>

#include <dune/grid/common/grid.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/refinement.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    template< class Grid, class RestrictProlongOperator >
    class AdaptRestrictProlongHandler
    {
      static const int dimension = Grid::dimension;

      typedef typename Grid::template Codim< 0 >::Entity Entity;
      typedef Dune::MakeableInterfaceObject< Entity > EntityObject;
      typedef typename EntityObject::ImplementationType EntityImp;

      typedef Alberta::ElementInfo< dimension > ElementInfo;
      typedef Alberta::Patch< dimension > Patch;

      Grid &grid_;
      RestrictProlongOperator &rpOp_;
      EntityObject father_;

    public:
      AdaptRestrictProlongHandler ( Grid &grid, RestrictProlongOperator &rpOp )
        : grid_( grid ),
          rpOp_( rpOp ),
          father_( EntityImp( grid_ ) )
      {}

      void restrictLocal ( const Patch &patch, int i )
      {
        ElementInfo fatherInfo = patch.elementInfo( i, grid_.levelProvider() );
        father_.impl().setElement( fatherInfo, 0 );
        rpOp_.preCoarsening( (const Entity &)father_ );
      }

      void prolongLocal ( const Patch &patch, int i )
      {
        ElementInfo fatherInfo = patch.elementInfo( i, grid_.levelProvider() );
        father_.impl().setElement( fatherInfo, 0 );
        rpOp_.postRefinement( (const Entity &)father_ );
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
