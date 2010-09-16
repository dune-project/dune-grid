// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDDATAHANDLE_HH
#define DUNE_ALU3DGRIDDATAHANDLE_HH

//- system includes
#include <iostream>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh>

//- local includes
#include "alu3dinclude.hh"

using std::endl;
using std::cout;
using std::flush;

namespace ALUGridSpace
{

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType, int codim >
  class GatherScatterBaseImpl
    : public GatherScatter
  {
  protected:
    const GridType & grid_;
    typedef typename GridType::template Codim<codim>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<codim>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType ImplElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    EntityType  & entity_;
    RealEntityType & realEntity_;

    DataCollectorType & dc_;

    const bool variableSize_;

    typedef typename GatherScatter :: ObjectStreamType ObjectStreamType;

    typedef typename DataCollectorType:: DataType DataType;
  public:
    //! Constructor
    GatherScatterBaseImpl(const GridType & grid, MakeableEntityType & en,
                          RealEntityType & realEntity , DataCollectorType & dc)
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

      // make sure partition type is set correct
      assert( elem.isGhost() == (entity_.partitionType() == Dune :: GhostEntity) );

      size_t size = getSize(str, entity_);
      // use normal scatter method
      dc_.scatter(str,entity_, size );
    }

    //! write Data of one element to stream
    void sendData ( ObjectStreamType & str , HElementType & elem )
    {
      // make sure element is contained in communication interface
      //assert( this->containsItem( elem ) );
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
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<0>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType ImplElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    typedef typename ImplTraits::template Codim< 1 >::InterfaceType HFaceType;

    typedef typename ImplTraits::template Codim< codim >::GhostInterfaceType HGhostType;
    typedef typename ImplTraits::template Codim< codim >::GhostImplementationType ImplGhostType;

    typedef typename ImplTraits::PllElementType PllElementType;

    EntityType& entity_;
    RealEntityType & realEntity_;

    // data handle
    DataCollectorType & dc_;

    const bool variableSize_;

    // used MessageBuffer
    typedef typename GatherScatter :: ObjectStreamType ObjectStreamType;

  public:
    //! Constructor
    GatherScatterBaseImpl(const GridType & grid, MakeableEntityType & en,
                          RealEntityType & realEntity , DataCollectorType & dc)
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
  template< class GridType, class DataCollectorType, int codim >
  class GatherScatterLeafData
    : public GatherScatterBaseImpl< GridType, DataCollectorType, codim >
  {
    enum { dim = GridType :: dimension };

    typedef GatherScatterBaseImpl<GridType,DataCollectorType,codim> BaseType;
    typedef typename GridType::template Codim<codim>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<codim>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType IMPLElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    typedef typename ImplTraits::template Codim< 1 >::InterfaceType HFaceType;

    typedef typename ImplTraits::template Codim< 0 >::GhostInterfaceType HGhostType;
    typedef typename ImplTraits::template Codim< 0 >::GhostImplementationType ImplGhostType;

    typedef typename ImplTraits::PllElementType PllElementType;

  public:
    //! Constructor
    GatherScatterLeafData(const GridType & grid, MakeableEntityType & en,
                          RealEntityType & realEntity , DataCollectorType & dc)
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
    typedef typename GridType::template Codim<codim>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<codim>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType IMPLElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    typedef typename ImplTraits::template Codim< 1 >::InterfaceType HFaceType;

    typedef typename ImplTraits::template Codim< 0 >::GhostInterfaceType HGhostType;
    typedef typename ImplTraits::template Codim< 0 >::GhostImplementationType ImplGhostType;

    typedef typename ImplTraits::PllElementType PllElementType;

    typedef typename GridType::LevelIndexSetImp LevelIndexSetImp;

    const LevelIndexSetImp & levelSet_;
    const int level_;
  public:
    //! Constructor
    GatherScatterLevelData(const GridType & grid, MakeableEntityType & en,
                           RealEntityType & realEntity , DataCollectorType & dc,
                           const LevelIndexSetImp & levelSet, const int level)
      : BaseType(grid,en,realEntity,dc) , levelSet_(levelSet) , level_(level)
    {}

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
    typedef typename GridType::template Codim<codim>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<codim>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType IMPLElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    typedef typename ImplTraits::template Codim< 1 >::InterfaceType HFaceType;

    typedef typename ImplTraits::template Codim< 0 >::GhostInterfaceType HGhostType;
    typedef typename ImplTraits::template Codim< 0 >::GhostImplementationType ImplGhostType;

    typedef typename ImplTraits::PllElementType PllElementType;

    typedef typename GridType::LevelIndexSetImp LevelIndexSetImp;

    const LevelIndexSetImp & levelSet_;
    const int level_;
  public:
    //! Constructor
    GatherScatterLevelData(const GridType & grid, MakeableEntityType & en,
                           RealEntityType & realEntity , DataCollectorType & dc,
                           const LevelIndexSetImp & levelSet, const int level)
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
      // if face level is not level_ then interior cannot be contained
      if(face.level() != level_) return false;

      typedef Gitter::helement_STI HElementType;
      typedef Gitter::hbndseg_STI HBndSegType;

      // check interior element here, might have a coarser level
      pair< HElementType *, HBndSegType * > p( (HElementType *)0, (HBndSegType *)0 );
      pll.getAttachedElement( p );
      assert( p.first );
      // check inside level
      bool contained = (p.first->level() == level_);
      assert( contained == this->containsItem( *p.first ));
      return contained;
    }

    // returns true, if ghost is contianed in set of comm interface
    bool containsGhost (const HFaceType & face, PllElementType & pll) const
    {
      // if face level is not level_ then ghost cannot be contained
      if(face.level() != level_) return false;
      // otherwise check ghost level
      return (pll.ghostLevel() == level_);
    }
  };
#endif // #if ALU3DGRID_PARALLEL

  //! the corresponding interface class is defined in bsinclude.hh
  template <class GridType, class DataCollectorType, class IndexOperatorType>
  class GatherScatterLoadBalance : public GatherScatter
  {
  protected:
    enum { codim = 0 };
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<0>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim< codim >::ImplementationType IMPLElementType;
    typedef typename ImplTraits::template Codim< codim >::InterfaceType HElementType;

    typedef typename ImplTraits::template Codim< 1 >::InterfaceType HFaceType;

    typedef typename ImplTraits::template Codim< 0 >::GhostInterfaceType HGhostType;
    typedef typename ImplTraits::template Codim< 0 >::GhostImplementationType ImplGhostType;

    typedef typename ImplTraits::PllElementType PllElementType;

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
    GatherScatterLoadBalance(GridType & grid, MakeableEntityType & en,
                             RealEntityType & realEntity , DataCollectorType & dc, IndexOperatorType & idxOp )
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
  template< class GridType, class AdaptDataHandle >
  class AdaptRestrictProlongImpl
    : public AdaptRestrictProlongType
  {
    GridType & grid_;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<0>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    //DofManagerType & dm_;
    AdaptDataHandle &rp_;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::HElementType HElementType;
    typedef typename ImplTraits::HBndSegType HBndSegType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

  public:
    //! Constructor
    AdaptRestrictProlongImpl ( GridType &grid,
                               MakeableEntityType &f, RealEntityType &rf,
                               MakeableEntityType &s, RealEntityType &rs,
                               AdaptDataHandle &rp )
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
      realFather_.setElement( elem );
      rp_.preCoarsening( reFather_ );

      // reset refinement marker
      elem.resetRefinedTag();
      return 0;
    }

    //! prolong data, elem is the father
    int postRefinement ( HElementType & elem )
    {
      realFather_.setElement( elem );
      rp_.postRefinement( reFather_ );

      // resert refinement markers
      elem.resetRefinedTag();
      for( HElementType *son = elem.down(); son != 0; son = son->next() )
        son->resetRefinedTag();

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



  template< class GridType, class AdaptDataHandle, class GlobalIdSetImp >
  class AdaptRestrictProlongGlSet
    : public AdaptRestrictProlongImpl< GridType, AdaptDataHandle >
  {
    typedef AdaptRestrictProlongImpl< GridType, AdaptDataHandle > BaseType;
    GlobalIdSetImp & set_;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<0>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::HElementType HElementType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

  public:
    //! Constructor
    AdaptRestrictProlongGlSet ( GridType &grid,
                                MakeableEntityType &f, RealEntityType &rf,
                                MakeableEntityType &s, RealEntityType &rs,
                                AdaptDataHandle &rp,
                                GlobalIdSetImp & set )
      : BaseType( grid, f, rf, s, rs, rp ),
        set_( set )
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
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef Dune :: MakeableInterfaceObject<
        typename GridType::template Codim<0>::Entity> MakeableEntityType;
    typedef typename MakeableEntityType :: ImplementationType RealEntityType;

    typedef typename GridType::Traits::LeafIndexSet LeafIndexSetType;

    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    DofManagerType & dm_;

    typedef typename GridType::MPICommunicatorType Comm;

    typedef Dune::ALU3dImplTraits< GridType::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::HElementType HElementType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

    int newMemSize_;
  public:
    //! Constructor
    LoadBalanceElementCount (GridType & grid,
                             MakeableEntityType & f, RealEntityType & rf, MakeableEntityType & s, RealEntityType & rs,DofManagerType & dm)
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

} // end namespace ALUGridSpace

#endif // #ifndef DUNE_ALU3DGRIDDATAHANDLE_HH
