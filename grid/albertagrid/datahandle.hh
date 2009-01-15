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
      EntityObject child_;
      ElementInfo fatherInfo_;
      ElementInfo childInfo_;
      RestrictProlongOperator &rpOp_;
      int maxlevel_;

    public:
      //! Constructor
      AdaptRestrictProlongHandler ( Grid &grid, RestrictProlongOperator &rpOp )
        : grid_( grid ),
          father_( EntityImp( grid_ ) ),
          child_( EntityImp( grid_ ) ),
          fatherInfo_( ElementInfo::createFake() ),
          childInfo_( ElementInfo::createFake() ),
          rpOp_( rpOp ),
          maxlevel_( -1 )
      {
        AlbertHelp::makeEmptyElInfo< dimension, dimension >( &(fatherInfo_.elInfo()) );
        AlbertHelp::makeEmptyElInfo< dimension, dimension >( &(childInfo_.elInfo()) );

        fatherInfo_.elInfo().mesh = grid.getMesh();
        childInfo_.elInfo().mesh = grid.getMesh();
      }

      //! restrict data , elem is always the father
      void preCoarsening ( Element *element )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        const int level = grid_.getLevelOfElement( element );

        fatherInfo_.elInfo().el = element;
        fatherInfo_.elInfo().level = level;
        Grid::getRealImplementation( father_ ).setElement( fatherInfo_, 0 );

#if DUNE_ALBERTA_VERSION >= 0x201
        childInfo_.elInfo().parent = &(fatherInfo_.elInfo());
#else
        childInfo_.elInfo().parent = element;
#endif
        childInfo_.elInfo().level = level+1;

        for( int i = 0; i < 2; ++i )
        {
          childInfo_.elInfo().el = element->child[ i ];
          assert( childInfo_.elInfo().el != NULL );
          Grid::getRealImplementation( child_ ).setElement( childInfo_, 0 );
          rpOp_.restrictLocal( father_, child_, (i == 0) );
        }

        maxlevel_ = std::max( level, maxlevel_ );
      }

      //! prolong data, elem is the father
      void postRefinement ( Alberta::Element *element )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        const int level = grid_.getLevelOfElement( element );

        fatherInfo_.elInfo().el = element;
        fatherInfo_.elInfo().level = level;

        Grid::getRealImplementation( father_ ).setElement( fatherInfo_, 0 );

#if DUNE_ALBERTA_VERSION >= 0x201
        childInfo_.elInfo().parent = &(fatherInfo_.elInfo());
#else
        childInfo_.elInfo().parent = element;
#endif
        childInfo_.elInfo().level = level+1;

        for( int i = 0; i < 2; ++i )
        {
          childInfo_.elInfo().el = element->child[ i ];
          assert( childInfo_.elInfo().el != NULL );
          Grid::getRealImplementation( child_ ).setElement( childInfo_, 0 );
          rpOp_.prolongLocal( father_, child_, (i == 0) );
        }

        maxlevel_ = std::max( level, maxlevel_ );
      }

#if 0
      int maxLevel () const
      {
        return maxlevel_;
      }
#endif
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
