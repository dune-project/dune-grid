// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDDATAHANDLE_HH
#define DUNE_ALBERTAGRIDDATAHANDLE_HH

//- system includes
#include <iostream>

#include <dune/grid/common/grid.hh>
#include <dune/grid/albertagrid/albertaheader.hh>

namespace AlbertHelp {

  /////////////////////////////////////////////////////////////////
  //
  //  --AdaptRestrictProlong
  //
  /////////////////////////////////////////////////////////////////
  template <class GridType , class RestrictProlongOperatorType >
  class AdaptRestrictProlongHandler
  {
    GridType & grid_;
    enum { dim = GridType :: dimension };
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    RestrictProlongOperatorType & rp_;

    ALBERTA EL_INFO vati_;
    ALBERTA EL_INFO sohn_;

    int maxlevel_;

  public:
    //! Constructor
    AdaptRestrictProlongHandler (GridType & grid,
                                 EntityType & f, RealEntityType & rf, EntityType & s, RealEntityType & rs
                                 , RestrictProlongOperatorType & rp)
      : grid_(grid)
        , reFather_(f)
        , reSon_(s)
        , realFather_(rf)
        , realSon_(rs)
        , rp_(rp)
        , vati_()
        , sohn_()
        , maxlevel_(-1)
    {
      makeEmptyElInfo<dim,dim>(&vati_);
      makeEmptyElInfo<dim,dim>(&sohn_);

      vati_.mesh = grid.getMesh();
      sohn_.mesh = grid.getMesh();
    }

    //! restrict data , elem is always the father
    void preCoarsening ( EL * elem )
    {
      // check pointers of mesh internal dof vecs
      // this call has to be done before all other things
      grid_.arrangeDofVec();

      vati_.el    = elem;
      int level   = grid_.getLevelOfElement( elem );
      vati_.level = level;
      ++level; // children level now

      realFather_.setElInfo(&vati_);
      sohn_.parent = elem;

      for(int ch=0; ch<2; ++ch)
      {
        EL * son = elem->child[ch];
        assert( son );

        sohn_.el    = son;
        sohn_.level = level;

        realSon_.setElInfo(&sohn_);

        rp_.restrictLocal(reFather_,reSon_,(ch == 0));
      }

      if(level > maxlevel_) maxlevel_ = level;
    }

    //! prolong data, elem is the father
    void postRefinement ( EL * elem )
    {
      // check pointers of mesh internal dof vecs
      // this call has to be done before all other things
      grid_.arrangeDofVec();

      vati_.el    = elem;
      int level   = grid_.getLevelOfElement( elem );
      vati_.level = level;
      ++level;
      realFather_.setElInfo(&vati_);
      sohn_.parent = elem;

      for(int ch=0; ch<2; ++ch)
      {
        EL * son = elem->child[ch];
        assert( son );

        sohn_.el    = son;
        sohn_.level = level;

        realSon_.setElInfo(&sohn_);

        rp_.prolongLocal(reFather_,reSon_,(ch == 0));
      }

      if(level > maxlevel_) maxlevel_ = level;
    }

    int maxLevel () const { return maxlevel_; }
  };

} // end namespace
#endif
