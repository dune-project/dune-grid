// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDDATAHANDLE_HH
#define DUNE_ALBERTAGRIDDATAHANDLE_HH

#if HAVE_ALBERTA

#include <iostream>

#include <dune/grid/common/grid.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/albertaextra.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

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

      Grid &grid_;
      EntityObject father_;
      RestrictProlongOperator &rpOp_;

    public:
      //! Constructor
      AdaptRestrictProlongHandler ( Grid &grid, RestrictProlongOperator &rpOp )
        : grid_( grid ),
          father_( EntityImp( grid_ ) ),
          rpOp_( rpOp )
      {}

      //! restrict data , elem is always the father
      void restrictLocal ( Element *father )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        const int level = grid_.getLevelOfElement( father );
        ElementInfo fatherInfo
          = ElementInfo::createFake( grid_.meshPointer(), father, level );
        Grid::getRealImplementation( father_ ).setElement( fatherInfo, 0 );

        rpOp_.restrictLocal( (const Entity &)father_ );
      }

      //! prolong data
      void prolongLocal ( Alberta::Element *father )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        const int level = grid_.getLevelOfElement( father );
        ElementInfo fatherInfo
          = ElementInfo::createFake( grid_.meshPointer(), father, level );
        Grid::getRealImplementation( father_ ).setElement( fatherInfo, 0 );

        rpOp_.prolongLocal( (const Entity &)father_ );
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
