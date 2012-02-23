// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDDATAHANDLE_HH
#define DUNE_ALU2DGRIDDATAHANDLE_HH

//- system includes
#include <iostream>

#include <dune/grid/common/adaptcallback.hh>

//- local includes
#include "alu2dinclude.hh"

using std::endl;
using std::cout;
using std::flush;

namespace ALU2DSPACENAME
{

  /////////////////////////////////////////////////////////////////
  //
  //  --AdaptRestrictProlong
  //
  /////////////////////////////////////////////////////////////////
  template< class GridType, class AdaptDataHandle >
  class AdaptRestrictProlong2dImpl
    : public AdaptRestrictProlong2d ALU2DDIMWORLD( GridType::dimensionworld, GridType::elementType )
  {
    GridType & grid_;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;
    typedef typename Dune::ALU2dImplTraits< GridType::dimensionworld, GridType::elementType >::HElementType HElementType ;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    AdaptDataHandle &rp_;

    int maxlevel_;


  public:
    //! Constructor
    AdaptRestrictProlong2dImpl ( GridType &grid,
                                 EntityType &f, RealEntityType &rf,
                                 EntityType &s, RealEntityType &rs,
                                 AdaptDataHandle &rp )
      : grid_(grid)
        , reFather_(f)
        , reSon_(s)
        , realFather_(rf)
        , realSon_(rs)
        , rp_(rp)
        , maxlevel_(-1)
    {}

    virtual ~AdaptRestrictProlong2dImpl ()
    {}

    //! restrict data , elem is always the father
    int preCoarsening ( HElementType &father )
    {
      maxlevel_ = std::max( maxlevel_, father.level() );
      //father.resetRefinedTag();
      realFather_.setElement( father );
      rp_.preCoarsening( reFather_ );

      return 0;
    }

    //! prolong data, elem is the father
    int postRefinement ( HElementType &father )
    {
      maxlevel_ = std::max( maxlevel_, father.level()+1 );
      //father.resetRefinedTag();
      realFather_.setElement( father );
      rp_.postRefinement( reFather_ );

      return 0;
    }

    int maxLevel () const { return maxlevel_; }
  };

} // namespace ALU2DSPACENAME

#endif // #ifndef DUNE_ALU2DGRIDDATAHANDLE_HH
