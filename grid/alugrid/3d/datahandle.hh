// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDDATAHANDLE_HH
#define DUNE_ALU3DGRIDDATAHANDLE_HH

#include <iostream>
#include "alu3dinclude.hh"

using std::endl;
using std::cout;
using std::flush;

namespace ALUGridSpace {

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType , int codim >
  class GatherScatterBaseImpl : public GatherScatter
  {
  protected:
    const GridType & grid_;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<codim>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::ImplementationType ImplElementType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::InterfaceType HElementType;

    EntityType     & entity_;
    RealEntityType & realEntity_;

    DataCollectorType & dc_;

    const bool variableSize_;

    typedef typename GatherScatter :: ObjectStreamType ObjectStreamType;

    typename DataCollectorType::DataType tmp_;
  public:
    //! Constructor
    GatherScatterBaseImpl(const GridType & grid, EntityType & en, RealEntityType & realEntity , DataCollectorType & dc)
      : grid_(grid), entity_(en), realEntity_(realEntity) , dc_(dc)
        , variableSize_( ! dc_.fixedsize(EntityType::dimension,codim) )
    {}

    //! returns contains of dc_
    bool contains(int dim, int cd) const { return dc_.contains(dim,cd); }

    // returns true, if element is contained in set of comm interface
    // this method must be overlaoded by the impl classes
    virtual bool containsItem (const HElementType & elem) const = 0;

    //! this method is called from the dunePackAll method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is written to the ObjectStream
    void inlineData ( ObjectStreamType & str , HElementType & elem )
    {
      // set element and then start
      realEntity_.setElement(elem);
      //dc_.inlineData(str,entity_);
      //std::cout << "inline entity " <<  grid_.hierarchicIndexSet().index(entity_) << "\n";
      dc_.gather(str,entity_);
    }

    //! this method is called from the duneUnpackSelf method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is read from the ObjectStream
    void xtractData ( ObjectStreamType & str , HElementType & elem )
    {
      // set element and then start
      grid_.updateStatus();
      realEntity_.setElement(elem);
      /*
         dc_.xtractData(str,entity_);
       */
      //std::cout << "xtract entity " <<  grid_.hierarchicIndexSet().index(entity_) << "\n";
      dc_.scatter(str,entity_,dc_.size(entity_));
    }

    void setData ( ObjectStreamType & str , HElementType & elem )
    {
      // one of this should be either true
      assert( this->containsItem( elem ) || elem.isGhost() );

      // set element and then start
      realEntity_.setElement(elem);
      // contiansItem might fail for ghost elements, which is not a bug
      // containsItem just returns right values for interior elements
      assert( entity_.partitionType() == Dune :: GhostEntity );

      size_t size = getSize(str, entity_);
      // use normal scatter method
      dc_.scatter(str,entity_, size );
    }

    //! write Data of one element to stream
    void sendData ( ObjectStreamType & str , HElementType & elem )
    {
      assert( this->containsItem( elem ) );
      realEntity_.setElement(elem);

      // if varaible size, also send size
      if( variableSize_ )
      {
        size_t size = dc_.size( entity_ );
        str.write( size );
      }

      dc_.gather(str, entity_ );
    }

    //! read Data of one element from stream
    void recvData ( ObjectStreamType & str , HElementType & elem )
    {
      assert( this->containsItem( elem ) );
      realEntity_.setElement( elem );

      size_t size = getSize(str, entity_);
      dc_.scatter(str,entity_, size );
    }

  protected:
    size_t getSize(ObjectStreamType & str, EntityType & en)
    {
      if(variableSize_)
      {
        size_t size;
        str.read(size);
        return size;
      }
      else
        return dc_.size(en);
    }
  };

  //***********************************************************
  //
  //  --specialisation for codim 0
  //
  //***********************************************************

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType >
  class GatherScatterBaseImpl<GridType,DataCollectorType,0> : public GatherScatter
  {
  protected:
    enum { codim = 0 };
    const GridType & grid_;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::ImplementationType IMPLElementType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::InterfaceType HElementType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<1>::InterfaceType HFaceType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::GhostInterfaceType HGhostType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::GhostImplementationType ImplGhostType;

#if ALU3DGRID_PARALLEL
    typedef ALU3DSPACE ElementPllXIF_t PllElementType;
#else
    typedef HElementType PllElementType;
#endif

    EntityType     & entity_;
    RealEntityType & realEntity_;

    // data handle
    DataCollectorType & dc_;

    const bool variableSize_;

    // used MessageBuffer
    typedef typename GatherScatter :: ObjectStreamType ObjectStreamType;

  public:
    //! Constructor
    GatherScatterBaseImpl(const GridType & grid, EntityType & en, RealEntityType & realEntity , DataCollectorType & dc)
      : grid_(grid), entity_(en), realEntity_(realEntity)
        , dc_(dc) , variableSize_ ( ! dc_.fixedsize( EntityType :: dimension, codim ))
    {}

    // return true if dim,codim combination is contained in data set
    bool contains(int dim, int codim) const
    {
      return dc_.contains(dim,codim);
    }

    // return true if item might from entity belonging to data set
    virtual bool containsItem (const HElementType & elem) const
    {
      return elem.isLeafEntity();
    }

    // return true if item might from entity belonging to data set
    virtual bool containsItem (const HGhostType & ghost) const = 0;

    //! this method is called from the dunePackAll method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is written to the ObjectStream
    void inlineData ( ObjectStreamType & str , HElementType & elem )
    {
      // set element and then start
      realEntity_.setElement(elem);
      dc_.gather(str,entity_);
      //dc_.inlineData(str,entity_);
    }

    //! this method is called from the duneUnpackSelf method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is read from the ObjectStream
    void xtractData ( ObjectStreamType & str , HElementType & elem )
    {
      // set element and then start
      const_cast<GridType &> (grid_).updateStatus();
      realEntity_.setElement(elem);
      dc_.scatter(str,entity_,dc_.size(entity_));
      // dc_.xtractData(str,entity_);
    }

    //! write Data of one element to stream
    void sendData ( ObjectStreamType & str , const HElementType & elem )
    {
      assert( this->containsItem(elem) );
      realEntity_.setElement( const_cast<HElementType &> (elem) );

      if( variableSize_ )
      {
        size_t size = dc_.size( entity_ );
        str.write( size );
      }
      dc_.gather(str, entity_);
    }

    //! read Data of one element from stream
    void recvData ( ObjectStreamType & str , HGhostType & ghost )
    {
      assert( this->containsItem( ghost ) );
      // don't check containsElement here, because this

      // set ghost as entity
      ImplGhostType & gh = static_cast <ImplGhostType &> (ghost);
      realEntity_.setGhost( gh );

      size_t size = getSize(str , entity_ );
      dc_.scatter(str, entity_, size );
    }

    //! read Data of one element from stream
    void recvData ( ObjectStreamType & str , HElementType & elem )
    {
      assert( this->containsItem( elem ) );
      realEntity_.setElement( elem );

      size_t size = getSize(str, entity_);
      dc_.scatter(str, entity_, size);
    }
  protected:
    size_t getSize(ObjectStreamType & str, EntityType & en)
    {
      if(variableSize_)
      {
        size_t size;
        str.read(size);
        return size;
      }
      else
        return dc_.size(en);
    }
  };

#if ALU3DGRID_PARALLEL
  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType , int codim >
  class GatherScatterLeafData
    : public GatherScatterBaseImpl<GridType,DataCollectorType,codim>
  {
    typedef GatherScatterBaseImpl<GridType,DataCollectorType,codim> BaseType;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<codim>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::ImplementationType IMPLElementType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::InterfaceType HElementType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<1>::InterfaceType HFaceType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostInterfaceType HGhostType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostImplementationType ImplGhostType;

    typedef ALU3DSPACE ElementPllXIF_t PllElementType;

  public:
    //! Constructor
    GatherScatterLeafData(const GridType & grid, EntityType & en, RealEntityType & realEntity , DataCollectorType & dc)
      : BaseType(grid,en,realEntity,dc) {}

    // returns true, if element is contained in set of comm interface
    bool containsItem (const HElementType & elem) const
    {
      return elem.isLeafEntity();
    }

    // returns true, if element is contained in set of comm interface
    bool containsItem (const HGhostType & ghost) const
    {
      return ghost.isLeafEntity();
    }

    // returns true, if interior element is contained in set of comm interface
    bool containsInterior (const HFaceType & face, PllElementType & pll) const
    {
      return face.isInteriorLeaf();
    }

    // returns true, if ghost is contianed in set of comm interface
    bool containsGhost (const HFaceType & face , PllElementType & pll) const
    {
      return pll.ghostLeaf();
    }
  };

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType , int codim >
  class GatherScatterLevelData
    : public GatherScatterBaseImpl<GridType,DataCollectorType,codim>
  {
    typedef GatherScatterBaseImpl<GridType,DataCollectorType,codim> BaseType;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<codim>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::ImplementationType IMPLElementType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::InterfaceType HElementType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<1>::InterfaceType HFaceType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostInterfaceType HGhostType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostImplementationType ImplGhostType;

    typedef ALU3DSPACE ElementPllXIF_t PllElementType;
    typedef typename GridType::LevelIndexSetImp LevelIndexSetImp;

    const LevelIndexSetImp & levelSet_;
    const int level_;
  public:
    //! Constructor
    GatherScatterLevelData(const GridType & grid, EntityType & en, RealEntityType & realEntity ,
                           DataCollectorType & dc, const LevelIndexSetImp & levelSet, const int level)
      : BaseType(grid,en,realEntity,dc) , levelSet_(levelSet) , level_(level) {}

    // returns true, if element is contained in set of comm interface
    bool containsItem (const HElementType & elem) const
    {
      return levelSet_.containsIndex(codim, elem.getIndex() );
    }
  };


  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType>
  class GatherScatterLevelData<GridType,DataCollectorType,0>
    : public GatherScatterBaseImpl<GridType,DataCollectorType,0>
  {
    enum { codim = 0 };
    typedef GatherScatterBaseImpl<GridType,DataCollectorType,codim> BaseType;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<codim>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::ImplementationType IMPLElementType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<codim>::InterfaceType HElementType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<1>::InterfaceType HFaceType;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostInterfaceType HGhostType;
    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::
    template Codim<0>::GhostImplementationType ImplGhostType;

    typedef ALU3DSPACE ElementPllXIF_t PllElementType;
    typedef typename GridType::LevelIndexSetImp LevelIndexSetImp;

    const LevelIndexSetImp & levelSet_;
    const int level_;
  public:
    //! Constructor
    GatherScatterLevelData(const GridType & grid, EntityType & en, RealEntityType & realEntity ,
                           DataCollectorType & dc, const LevelIndexSetImp & levelSet, const int level)
      : BaseType(grid,en,realEntity,dc) , levelSet_(levelSet) , level_(level) {}

    // returns true, if element is contained in set of comm interface
    bool containsItem (const HElementType & elem) const
    {
      return levelSet_.containsIndex(codim, elem.getIndex() );
    }

    // returns true, if element is contained in set of comm interface
    bool containsItem (const HGhostType & ghost) const
    {
      assert( ghost.getGhost().first );
      return containsItem( * (ghost.getGhost().first) );
    }

    // returns true, if interior element is contained in set of comm interface
    bool containsInterior (const HFaceType & face, PllElementType & pll) const
    {
      if(face.level() != level_) return false;

      // check interior element here, might have a coarser level
      pair < Gitter::helement_STI * , Gitter::hbndseg_STI * > p (0,0);
      pll.getAttachedElement( p );
      assert( p.first );
      bool contained = (p.first->level() == level_);
      assert( contained == this->containsItem( *p.first ));
      return contained;
    }

    // returns true, if ghost is contianed in set of comm interface
    bool containsGhost (const HFaceType & face, PllElementType & pll) const
    {
      return (pll.ghostLevel() == level_);
    }
  };
#endif

  /////////////////////////////////////////////////////////////////
  //
  //  --AdaptRestrictProlong
  //
  /////////////////////////////////////////////////////////////////
  template <class GridType , class RestrictProlongOperatorType >
  class AdaptRestrictProlongImpl : public AdaptRestrictProlongType
  {
    GridType & grid_;
    typedef typename GridType::template Codim<0>::Entity Entity;
    typedef typename Entity :: ImplementationType EntityImp;
    typedef typename GridType::Traits::LeafIndexSet LeafIndexSetType;

    Entity & reFather_;
    Entity & reSon_;
    EntityImp & realFather_;
    EntityImp & realSon_;

    //DofManagerType & dm_;
    RestrictProlongOperatorType & rp_;

    int maxlevel_;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::PLLBndFaceType PLLBndFaceType;

  public:
    //! Constructor
    AdaptRestrictProlongImpl (GridType & grid,
                              Entity & f, EntityImp & rf, Entity & s, EntityImp & rs
                              , RestrictProlongOperatorType & rp)
      : grid_(grid)
        , reFather_(f)
        , reSon_(s)
        , realFather_(rf)
        , realSon_(rs)
        , rp_(rp)
        , maxlevel_(-1)
    {}

    virtual ~AdaptRestrictProlongImpl ()
    {}

    //! restrict data , elem is always the father
    int preCoarsening ( HElementType & elem )
    {
      // set element and then start
      HElementType * son = elem.down();
      if(elem.level() > maxlevel_) maxlevel_ = elem.level();

      //elem.resetRefinedTag();
      assert( son );

      realSon_.setElement(*son);
      realFather_.setElement(elem);
      rp_.restrictLocal(reFather_,reSon_,true);

      son = son->next();
      while( son )
      {
        realSon_.setElement(*son);
        rp_.restrictLocal(reFather_,reSon_,false);
        son = son->next();
      }
      return 0;
    }

    //! prolong data, elem is the father
    int postRefinement ( HElementType & elem )
    {
      // set element and then start
      HElementType * son = elem.down();
      assert( son );
      //elem.resetRefinedTag();

      realFather_.setElement(elem);
      realSon_.setElement(*son);
      int level = realSon_.level();
      if(level > maxlevel_) maxlevel_ = level;

      rp_.prolongLocal(reFather_,reSon_, false);

      son = son->next();
      while( son )
      {
        assert( son );

        realSon_.setElement(*son);
        rp_.prolongLocal(reFather_,reSon_, false);

        son = son->next();
      }
      return 0;
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    int preCoarsening ( HBndSegType & el )
    {
      /*
         // hier nur die Indices einfügen
         // da wir nur 4 ghost kenn macht das prolong keinen Sinn

         PLLBndFaceType & elem = static_cast<PLLBndFaceType &> (el);

         realFather_.setGhost(elem);
         dm_.insertNewIndex( reFather_ );

         // set element and then start
         PLLBndFaceType * son = static_cast<PLLBndFaceType *> (elem.down());
         while( son )
         {
         realSon_.setGhost(*son);
         dm_.removeOldIndex( reSon_ );
         son = static_cast<PLLBndFaceType *> (son->next());
         }
       */
      return 0;
    }


    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el )
    {
      /*
         // for ghost only indices are inserted, because we only know 4 of 8
         // elements and therefore no prolongation can be done
         PLLBndFaceType & elem = static_cast<PLLBndFaceType &> (el);

         realFather_.setGhost(elem);
         dm_.insertNewIndex( reFather_ );

         // set element and then start
         PLLBndFaceType * son = static_cast<PLLBndFaceType *> (elem.down());
         while( son )
         {
         assert( son );
         realSon_.setGhost(*son);
         dm_.insertNewIndex( reSon_ );
         son = static_cast<PLLBndFaceType *> (son->next());
         }
       */
      return 0;
    }
    int maxLevel () const { return maxlevel_; }
  };

  template <class GridType , class RestrictProlongOperatorType ,
      class GlobalIdSetImp >
  class AdaptRestrictProlongGlSet
    : public AdaptRestrictProlongImpl<GridType,RestrictProlongOperatorType>
  {
    typedef AdaptRestrictProlongImpl<GridType,RestrictProlongOperatorType> BaseType;
    GlobalIdSetImp & set_;
    typedef typename GridType::template Codim<0>::Entity Entity;
    typedef typename Entity :: ImplementationType EntityImp;

  public:
    //! Constructor
    AdaptRestrictProlongGlSet(GridType & grid,
                              Entity & f, EntityImp & rf, Entity & s, EntityImp & rs
                              , RestrictProlongOperatorType & rp
                              , GlobalIdSetImp & set )
      : BaseType(grid,f,rf,s,rs,rp)
        , set_(set)
    {}

    virtual ~AdaptRestrictProlongGlSet () {}

    //! restrict data , elem is always the father
    int preCoarsening ( HElementType & elem )
    {
      return BaseType :: preCoarsening (elem );
    }

    //! prolong data, elem is the father
    int postRefinement ( HElementType & elem )
    {
      set_.postRefinement( elem );
      return BaseType :: postRefinement(elem );
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    int preCoarsening ( HBndSegType & el )
    {
      return 0;
    }


    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el )
    {
      return 0;
    }

  };

  //***************************************************************
  //
  //  for load balancing
  //
  //***************************************************************

  template <class GridType, class DofManagerType>
  class LoadBalanceRestrictProlongImpl : public AdaptRestrictProlongType
  {
    GridType & grid_;
    typedef typename GridType:: EntityObject Entity;
    typedef typename Entity :: ImplementationType EntityImp;

    Entity & reFather_;
    EntityImp & realFather_;

    Entity & reSon_;
    EntityImp & realSon_;

    DofManagerType & dm_;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::PLLBndFaceType PLLBndFaceType;

    int newMemSize_;

  public:
    //! Constructor
    LoadBalanceRestrictProlongImpl (GridType & grid, Entity & f, EntityImp &
                                    rf, Entity & s, EntityImp & rs, DofManagerType & dm )
      : grid_(grid), reFather_(f), realFather_(rf) , reSon_(s), realSon_(rs)
        , dm_(dm), newMemSize_ (0) {}

    virtual ~LoadBalanceRestrictProlongImpl () {};

    //! restrict data , elem is always the father
    int postRefinement ( HElementType & elem )
    {
      realFather_.setElement(elem);
      dm_.removeOldIndex( reFather_ );

      HElementType * son = elem.down();
      while( son )
      {
        realSon_.setElement(*son);
        dm_.insertNewIndex( reSon_ );
        son = son->next();
        newMemSize_ ++ ;
      }
      return 0;
    }

    //! prolong data, elem is the father
    int preCoarsening ( HElementType & elem )
    {
      realFather_.setElement(elem);
      dm_.insertNewIndex( reFather_ );

      HElementType * son = elem.down();
      while( son )
      {
        realSon_.setElement(*son);
        dm_.removeOldIndex( reSon_ );
        son = son->next();
      }
      return 0;
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    int preCoarsening ( HBndSegType & el )
    {
      PLLBndFaceType & elem = static_cast<PLLBndFaceType &> (el);

      realFather_.setGhost(elem);
      dm_.insertNewIndex( reFather_ );

      // set element and then start
      PLLBndFaceType * son = static_cast<PLLBndFaceType *> (elem.down());
      while( son )
      {
        realSon_.setGhost(*son);
        dm_.removeOldIndex( reSon_ );
        son = static_cast<PLLBndFaceType *> (son->next());
      }
      return 0;
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    int postRefinement ( HBndSegType & el )
    {
      PLLBndFaceType & elem = static_cast<PLLBndFaceType &> (el);

      PLLBndFaceType * vati = static_cast<PLLBndFaceType *> ( elem.up());
      if( vati )
      {
        realFather_.setGhost(*vati);
        dm_.removeOldIndex( reFather_ );
      }

      realFather_.setGhost(elem);
      dm_.insertNewIndex( reFather_ );
      newMemSize_ ++;

      // set element and then start
      PLLBndFaceType * son = static_cast<PLLBndFaceType *> (elem.down());
      while( son )
      {
        realSon_.setGhost(*son);
        dm_.insertNewIndex( reSon_ );
        son = static_cast<PLLBndFaceType *> (son->next());
        newMemSize_ ++;
      }
      return 0;
    }

    int newElements () const { return newMemSize_; }
  };

} // end namespace

#endif
