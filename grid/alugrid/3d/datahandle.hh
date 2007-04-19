// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDDATAHANDLE_HH
#define DUNE_ALU3DGRIDDATAHANDLE_HH

//- system includes
#include <iostream>

//- local includes
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

    // set elem to realEntity
    virtual void setElement(const HElementType & elem) = 0;

    void setData ( ObjectStreamType & str , HElementType & elem )
    {
      // one of this should be either true
      assert( this->containsItem( elem ) || elem.isGhost() );

      // set element and then start
      setElement(elem);

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
      setElement(elem);

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
      setElement( elem );

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
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
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

    //! write Data of one element to stream
    void sendData ( ObjectStreamType & str , const HElementType & elem )
    {
      assert( this->containsItem(elem) );
      realEntity_.setElement( const_cast<HElementType &> (elem) );

      // write size in case of variable size
      writeSize( str, entity_);
      // gather data
      dc_.gather(str, entity_);
    }

    //! write Data of one ghost element to stream
    void sendData ( ObjectStreamType & str , const HGhostType& ghost)
    {
      assert( this->containsItem( ghost ) );

      // set ghost as entity
      realEntity_.setGhost( const_cast <HGhostType &> (ghost) );

      // write size in case of variable size
      writeSize( str, entity_);
      // gather data
      dc_.gather(str, entity_);
    }

    //! read Data of one element from stream
    void recvData ( ObjectStreamType & str , HElementType & elem )
    {
      assert( this->containsItem( elem ) );
      realEntity_.setElement( elem );

      size_t size = getSize(str, entity_);
      dc_.scatter(str, entity_, size);
    }

    //! read Data of one element from stream
    void recvData ( ObjectStreamType & str , HGhostType & ghost )
    {
      assert( this->containsItem( ghost ) );

      // set ghost as entity
      realEntity_.setGhost( ghost );

      size_t size = getSize(str , entity_ );
      dc_.scatter(str, entity_, size );
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

    // write variable size to stream
    void writeSize(ObjectStreamType & str, EntityType & en)
    {
      if( variableSize_ )
      {
        size_t size = dc_.size( en );
        str.write( size );
      }
    }
  };

#if ALU3DGRID_PARALLEL
  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType , int codim >
  class GatherScatterLeafData
    : public GatherScatterBaseImpl<GridType,DataCollectorType,codim>
  {
    enum { dim = GridType :: dimension };

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
      : BaseType(grid,en,realEntity,dc)
    {
      // if leaf vertices are communicated,
      // make sure that vertex list is up2date
      // but only do this, if vertex data contained,
      // because the list update is expensive
      if( (codim == 3) && dc.contains(dim,codim) )
      {
        // call of this method forces update of list,
        // if list is not up to date
        grid.getLeafVertexList();
      }
    }

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

    // set elem to realEntity
    void setElement(const HElementType & elem)
    {
      this->realEntity_.setElement(elem);
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

    // set elem to realEntity
    void setElement(const HElementType & elem)
    {
      this->realEntity_.setElement(elem,level_);
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

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType, class IndexOperatorType>
  class GatherScatterLoadBalance : public GatherScatter
  {
  protected:
    enum { codim = 0 };
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
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
    GridType & grid_;

    EntityType     & entity_;
    RealEntityType & realEntity_;

    // data handle
    DataCollectorType & dc_;
    IndexOperatorType & idxOp_;

    // used MessageBuffer
    typedef typename GatherScatter :: ObjectStreamType ObjectStreamType;

  public:
    //! Constructor
    GatherScatterLoadBalance(GridType & grid, EntityType & en, RealEntityType & realEntity , DataCollectorType & dc, IndexOperatorType & idxOp )
      : grid_(grid), entity_(en), realEntity_(realEntity)
        , dc_(dc) , idxOp_(idxOp)
    {}

    // return true if dim,codim combination is contained in data set
    bool contains(int dim, int codim) const
    {
      return true;
    }

    //! this method is called from the dunePackAll method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is written to the ObjectStream
    void inlineData ( ObjectStreamType & str , HElementType & elem )
    {
      str.write(grid_.maxLevel());
      // set element and then start
      assert( elem.level () == 0 );
      realEntity_.setElement(elem);
      dc_.inlineData(str,entity_);
    }

    //! this method is called from the duneUnpackSelf method of the corresponding
    //! Macro element class of the BSGrid, see gitter_dune_pll*.*
    //! here the data is read from the ObjectStream
    void xtractData ( ObjectStreamType & str , HElementType & elem )
    {
      assert( elem.level () == 0 );
      int mxl;
      str.read(mxl);
      // set element and then start
      grid_.setMaxLevel(mxl);

      // reserve memory for new elements
      size_t elChunk = idxOp_.newElements();
      assert( elChunk > 0 );

      realEntity_.setElement(elem);
      dc_.xtractData(str,entity_, elChunk);
    }

    //! call compress on data
    void compress ()
    {
      dc_.compress();
    }
  };

  /////////////////////////////////////////////////////////////////
  //
  //  --AdaptRestrictProlong
  //
  /////////////////////////////////////////////////////////////////
  template <class GridType , class RestrictProlongOperatorType >
  class AdaptRestrictProlongImpl : public AdaptRestrictProlongType
  {
    GridType & grid_;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    //DofManagerType & dm_;
    RestrictProlongOperatorType & rp_;

    typedef typename Dune::ALU3dImplTraits<GridType::elementType>::PLLBndFaceType PLLBndFaceType;

  public:
    //! Constructor
    AdaptRestrictProlongImpl (GridType & grid,
                              EntityType & f, RealEntityType & rf, EntityType & s, RealEntityType & rs
                              , RestrictProlongOperatorType & rp)
      : grid_(grid)
        , reFather_(f)
        , reSon_(s)
        , realFather_(rf)
        , realSon_(rs)
        , rp_(rp)
    {}

    virtual ~AdaptRestrictProlongImpl ()
    {}

    //! restrict data , elem is always the father
    int preCoarsening ( HElementType & elem )
    {
      // set element and then start
      HElementType * son = elem.down();

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

      // reset refinement marker
      elem.resetRefinedTag();
      return 0;
    }

    //! prolong data, elem is the father
    int postRefinement ( HElementType & elem )
    {
      // set element and then start
      HElementType * son = elem.down();
      assert( son );

      // reset refinement marker
      elem.resetRefinedTag();

      realFather_.setElement(elem);
      realSon_.setElement(*son);

      rp_.prolongLocal(reFather_,reSon_, false);

      son = son->next();
      while( son )
      {
        assert( son );

        realSon_.setElement(*son);
        rp_.prolongLocal(reFather_,reSon_, false);

        // reset refinement tag
        son->resetRefinedTag();

        son = son->next();
      }
      return 0;
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

  template <class GridType , class RestrictProlongOperatorType ,
      class GlobalIdSetImp >
  class AdaptRestrictProlongGlSet
    : public AdaptRestrictProlongImpl<GridType,RestrictProlongOperatorType>
  {
    typedef AdaptRestrictProlongImpl<GridType,RestrictProlongOperatorType> BaseType;
    GlobalIdSetImp & set_;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

  public:
    //! Constructor
    AdaptRestrictProlongGlSet(GridType & grid,
                              EntityType & f, RealEntityType & rf, EntityType & s, RealEntityType & rs
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

  };

  // this class is for counting the tree depth of the
  // element when unpacking data from load balance
  template <class GridType , class DofManagerType>
  class LoadBalanceElementCount : public AdaptRestrictProlongType
  {
    GridType & grid_;
    typedef Dune :: MakeableInterfaceObject<typename GridType::template Codim<0>::Entity> EntityType;
    typedef typename EntityType :: ImplementationType RealEntityType;

    typedef typename GridType::Traits::LeafIndexSet LeafIndexSetType;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    DofManagerType & dm_;
    typedef typename  Dune::ALU3dImplTraits<GridType::elementType>::PLLBndFaceType PLLBndFaceType;

    int newMemSize_;
  public:
    //! Constructor
    LoadBalanceElementCount (GridType & grid,
                             EntityType & f, RealEntityType & rf, EntityType & s, RealEntityType & rs,DofManagerType & dm)
      : grid_(grid)
        , reFather_(f)
        , reSon_(s)
        , realFather_(rf)
        , realSon_(rs)
        , dm_(dm)
        , newMemSize_ (1) // we have at least one element (the macro element)
    {}

    virtual ~LoadBalanceElementCount () {};

    //! restrict data , elem is always the father
    int postRefinement ( HElementType & elem )
    {
      // when called for a macro element, then a new tree is starting
      // set to 1 because for only macro elements this method is not called
      if( elem.level() == 0 ) newMemSize_ = 1;

      for( HElementType * son = elem.down() ; son ; son= son->next())
      {
        ++ newMemSize_;
      }
      return 0;
    }

    //! prolong data, elem is the father
    int preCoarsening ( HElementType & elem )
    {
      return 0;
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    int preCoarsening ( HBndSegType & el )
    {
      return 0;
    }

    //! restrict data , elem is always the father
    //! this method is for ghost elements
    //! we need the ghost method because data is only inlined for interior
    //! elements, but we have to arange the ghost indices
    int postRefinement ( HBndSegType & el )
    {
      return 0;
    }

    int newElements () const { return newMemSize_; }
  };

} // end namespace
#endif
