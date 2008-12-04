// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDDATAHANDLE_HH
#define DUNE_ALBERTAGRIDDATAHANDLE_HH

#if HAVE_ALBERTA

//- system includes
#include <iostream>

#include <dune/grid/common/grid.hh>
#include <dune/grid/albertagrid/albertaheader.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

namespace Dune
{

  namespace Alberta
  {

    template <class GridType , class RestrictProlongOperatorType >
    class AdaptRestrictProlongHandler
    {
      static const int dimension = GridType::dimension;

      GridType & grid_;

      typedef Dune::MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
      typedef typename EntityType :: ImplementationType RealEntityType;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

      EntityType & reFather_;
      EntityType & reSon_;
      RealEntityType & realFather_;
      RealEntityType & realSon_;

      RestrictProlongOperatorType &rp_;

      ElementInfo fatherInfo_;
      ElementInfo childInfo_;

      int maxlevel_;

    public:
      //! Constructor
      AdaptRestrictProlongHandler ( GridType &grid,
                                    EntityType &f, RealEntityType &rf,
                                    EntityType &s, RealEntityType &rs,
                                    RestrictProlongOperatorType &rp )
        : grid_( grid ),
          reFather_( f ),
          reSon_( s ),
          realFather_( rf ),
          realSon_( rs ),
          rp_( rp ),
          fatherInfo_( ElementInfo::createFake() ),
          childInfo_( ElementInfo::createFake() ),
          maxlevel_( -1 )
      {
        AlbertHelp::makeEmptyElInfo< dimension, dimension >( &(fatherInfo_.elInfo()) );
        AlbertHelp::makeEmptyElInfo< dimension, dimension >( &(childInfo_.elInfo()) );

        fatherInfo_.elInfo().mesh = grid.getMesh();
        childInfo_.elInfo().mesh = grid.getMesh();
      }

      //! restrict data , elem is always the father
      void preCoarsening ( EL * elem )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        int level   = grid_.getLevelOfElement( elem );

        fatherInfo_.elInfo().el = elem;
        fatherInfo_.elInfo().level = level;

        realFather_.setElement( fatherInfo_, 0 );

#if DUNE_ALBERTA_VERSION >= 0x201
        childInfo_.elInfo().parent = &(fatherInfo_.elInfo());
#else
        childInfo_.elInfo().parent = elem;
#endif
        childInfo_.elInfo().level = level+1;

        for( int i = 0; i < 2; ++i )
        {
          EL *child = elem->child[ i ];
          assert( child != NULL );

          childInfo_.elInfo().el = child;

          realSon_.setElement( childInfo_, 0 );

          rp_.restrictLocal( reFather_, reSon_, (i == 0) );
        }

        if( level+1 > maxlevel_ )
          maxlevel_ = level;
      }

      //! prolong data, elem is the father
      void postRefinement ( EL * elem )
      {
        // check pointers of mesh internal dof vecs
        // this call has to be done before all other things
        grid_.arrangeDofVec();

        int level   = grid_.getLevelOfElement( elem );

        fatherInfo_.elInfo().el = elem;
        fatherInfo_.elInfo().level = level;

        realFather_.setElement( fatherInfo_, 0 );

#if DUNE_ALBERTA_VERSION >= 0x201
        childInfo_.elInfo().parent = &(fatherInfo_.elInfo());
#else
        childInfo_.elInfo().parent = elem;
#endif
        childInfo_.elInfo().level = level+1;

        for( int i = 0; i < 2; ++i )
        {
          EL *child = elem->child[ i ];
          assert( child );

          childInfo_.elInfo().el = child;

          realSon_.setElement( childInfo_, 0 );

          rp_.prolongLocal( reFather_, reSon_, (i == 0) );
        }

        if( level+1 > maxlevel_ )
          maxlevel_ = level;
      }

      int maxLevel () const
      {
        return maxlevel_;
      }
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
