// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Dune includes
#include <dune/common/stdstreams.hh>

// Local includes
#include "entity.hh"
#include "iterator.hh"
#include "myautoptr.hh"
#include "datahandle.hh"

namespace Dune {


  template <class GridType >
  inline void ALU3dGridVertexList ::
  setupVxList(const GridType & grid, int level)
  {
    // iterates over grid elements of given level and adds all vertices to
    // given list
    enum { codim = 3 };

    VertexListType & vxList = vertexList_;

    unsigned int vxsize = grid.hierarchicIndexSet().size(codim);
    if( vxList.size() < vxsize ) vxList.reserve(vxsize);
    Array<int> visited_(vxsize);

    for(unsigned int i=0; i<vxsize; i++)
    {
      visited_[i] = 0;
    }

    vxList.resize(0);

    const ALU3dGridElementType elType = GridType:: elementType;

    typedef ALU3DSPACE ALU3dGridLevelIteratorWrapper<0,Dune::All_Partition>  ElementLevelIteratorType;
    typedef typename ElementLevelIteratorType :: val_t val_t;

    typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;
    typedef ALU3DSPACE VertexType VertexType;

    enum { nVx = ElementTopologyMapping < elType > :: numVertices };

    ElementLevelIteratorType it ( grid, *this , level, grid.nlinks() );

    int count = 0;
    for( it.first(); !it.done() ; it.next())
    {
      val_t & item = it.item();

      IMPLElementType * elem = 0;
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first);
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first);

      assert( elem );

      for(int i=0; i<nVx; ++i)
      {
        VertexType * vx = elem->myvertex(i);
        assert( vx );

        // insert only interior and border vertices
        if( vx->isGhost() ) continue;

        const int idx = vx->getIndex();
        if(visited_[idx] == 0)
        {
          vxList.push_back(vx);
          ++count;
        }
        visited_[idx] = 1;
      }
    }
    assert( count == (int) vxList.size());;
    up2Date_ = true;
  }


  //--Grid
  //template <int dim, int dimworld, ALU3dGridElementType elType>
  //const ALU3dGridElementType
  //ALU3dGrid<dim, dimworld, elType>::elementType = elType;

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::
  ALU3dGrid(const std::string macroTriangFilename
#if ALU3DGRID_PARALLEL
            , const MPI_Comm mpiComm
#endif
            )
    : mygrid_ (0)
#if ALU3DGRID_PARALLEL
      , mpAccess_(mpiComm)
      , myRank_( mpAccess_.myrank() )
      , ccobj_(mpiComm)
#else
      , myRank_(-1)
#endif
      , maxlevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_ (*this)
      , globalIdSet_(0), localIdSet_(*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , sizeCache_ (0)
      , ghostElements_(0)
  {
    makeGeomTypes();

    mygrid_ = new ALU3DSPACE GitterImplType (macroTriangFilename.c_str()
#if ALU3DGRID_PARALLEL
                                             , mpAccess_
#endif
                                             );
    assert(mygrid_ != 0);

#if ALU3DGRID_PARALLEL
    dverb << "************************************************\n";
    dverb << "Created grid on p=" << mpAccess_.myrank() << "\n";
    dverb << "************************************************\n";
#endif
    this->checkMacroGrid ();

    postAdapt();
    calcExtras();

    if (size(0)>0)
    {
      std::cout << "Created ALU3dGrid from macro grid file '"
                << macroTriangFilename << "'. \n\n";
    }
    /*
       else
       {
       std::cout << "Created empty ALU3dGrid. \n\n";
       }
     */
  }

#if ALU3DGRID_PARALLEL
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::ALU3dGrid(const MPI_Comm mpiComm)
    : mygrid_ (0)
      , mpAccess_(mpiComm)
      , myRank_( mpAccess_.myrank() )
      , ccobj_(mpiComm)
      , maxlevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_ (*this)
      , globalIdSet_(0), localIdSet_(*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , sizeCache_ (0)
      , ghostElements_(0)
  {
    makeGeomTypes();
  }
#else
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::ALU3dGrid(int myrank)
    : mygrid_ (0)
      , myRank_(myrank)
      , maxlevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_ (*this)
      , globalIdSet_ (0)
      , localIdSet_ (*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , ghostElements_(0)
  {
    makeGeomTypes();
  }
#endif

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::ALU3dGrid(const ALU3dGrid<dim, dimworld, elType> & g)
    : mygrid_ (0)
#if ALU3DGRID_PARALLEL
      , mpAccess_(g.mpAccess_)
#endif
      , myRank_(-1)
      , ccobj_(g.ccobj_)
      , maxlevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_(*this)
      , globalIdSet_ (0)
      , localIdSet_ (*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , sizeCache_ (0)
      , ghostElements_(0)
  {
    DUNE_THROW(GridError,"Do not use copy constructor of ALU3dGrid! \n");
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType> & ALU3dGrid<dim, dimworld, elType>::operator = (const ALU3dGrid<dim, dimworld, elType> & g)
  {
    DUNE_THROW(GridError,"Do not use assignment operator of ALU3dGrid! \n");
    return (*this);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::~ALU3dGrid()
  {
    for(unsigned int i=0; i<levelIndexVec_.size(); i++) delete levelIndexVec_[i];
    delete globalIdSet_; globalIdSet_ = 0;
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;
    delete mygrid_; mygrid_ = 0;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::size(int level, int codim) const
  {
    // if we dont have this level return 0
    if( (level > maxlevel_) || (level < 0) ) return 0;

    assert( codim >= 0);
    assert( codim < dim+1 );

    //assert( levelIndexSet(level).size(codim,this->geomTypes(codim)[0]) ==
    //   sizeCache_->size(level,codim) );
    assert( sizeCache_ );
    return sizeCache_->size(level,codim);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::makeGeomTypes()
  {
    if(elType == tetra)
    {
      // stored is the dim, where is the codim
      for(int i=dim; i>= 0; i--)
        geomTypes_[dim-i][0] = GeometryType(GeometryType::simplex,i);
      return;
    }
    if(elType == hexa)
    {
      // stored is the dim, where is the codim
      for(int i=dim; i>= 0; i--)
        geomTypes_[dim-i][0] = GeometryType(GeometryType::cube,i);
      return;
    }
    DUNE_THROW(GridError,"Geometrytype not implemented!");
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::size(int level, GeometryType type) const
  {
    if(elType == tetra && !type.isSimplex()) return 0;
    if(elType == hexa  && !type.isCube   ()) return 0;
    return size(level,dim-type.dim());
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::size(int codim) const
  {
    assert( codim >= 0 );
    assert( codim < dim +1 );

    assert( sizeCache_ );
    return sizeCache_->size(codim);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::size(GeometryType type) const
  {
    if(elType == tetra && !type.isSimplex()) return 0;
    if(elType == hexa  && !type.isCube   ()) return 0;
    return size(dim-type.dim());
  }

  // calc all necessary things that might have changed
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::updateStatus()
  {
    calcMaxlevel();
    calcExtras();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::calcMaxlevel()
  {
    maxlevel_ = 0;
    ALU3DSPACE BSLeafIteratorMaxLevel w (myGrid()) ;
    for (w->first () ; ! w->done () ; w->next ())
    {
      if(w->item().level() > maxlevel_ ) maxlevel_ = w->item().level();
    }
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::calcExtras()
  {
    if(sizeCache_) delete sizeCache_;
    bool isSimplex = (elType == tetra) ? true : false;
    sizeCache_ = new SizeCacheType (*this,isSimplex,!isSimplex,true);

    for(size_t i=0; i<levelIndexVec_.size(); i++)
    {
      if(levelIndexVec_[i]) (*(levelIndexVec_[i])).calcNewIndex();
    }

    if(leafIndexSet_) leafIndexSet_->calcNewIndex();

    // update id set, i.e. insert new elements
    //if(globalIdSet_) globalIdSet_->updateIdSet();
#ifndef NDEBUG
    //if(globalIdSet_) globalIdSet_->uniquenessCheck();
#endif

    for(size_t i=0; i<MAXL; i++)
    {
      vertexList_[i].unsetUp2Date();
      levelEdgeList_[i].unsetUp2Date();
    }
#if ALU3DGRID_PARALLEL
    for(int i=0; i<dim; ++i)
    {
      ghostLeafList_[i].unsetUp2Date();
      for(size_t l=0; l<MAXL; ++l) ghostLevelList_[i][l].unsetUp2Date();
    }
#endif
    coarsenMarked_ = 0;
    refineMarked_  = 0;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::global_size(int codim) const
  {
    // return actual size of hierarchical index set
    // this is always up to date
    // maxIndex is the largest index used + 1
    assert( mygrid_ );
    return (*mygrid_).indexManager(codim).getMaxIndex();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::hierSetSize(int codim) const
  {
    // return actual size of hierarchical index set
    assert( mygrid_ );
    return (*mygrid_).indexManager(codim).getMaxIndex();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline const typename ALU3dGrid<dim, dimworld, elType>::Traits :: LeafIndexSet &
  ALU3dGrid<dim, dimworld, elType>::leafIndexSet() const
  {
    if(!leafIndexSet_) leafIndexSet_ = new LeafIndexSetImp ( *this );
    return *leafIndexSet_;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline const typename ALU3dGrid<dim, dimworld, elType>::Traits :: LevelIndexSet &
  ALU3dGrid<dim, dimworld, elType>::levelIndexSet( int level ) const
  {
    // check if level fits in vector
    assert( level >= 0 );
    assert( level < (int) levelIndexVec_.size() );

    if( levelIndexVec_[level] == 0 )
      levelIndexVec_[level] = new LevelIndexSetImp ( *this , level );
    return *(levelIndexVec_[level]);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::maxLevel() const
  {
    return maxlevel_;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3DSPACE GitterImplType & ALU3dGrid<dim, dimworld, elType>::myGrid() const
  {
    assert( mygrid_ );
    return *mygrid_;
  }

  // lbegin methods
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU3dGrid<dim, dimworld, elType>::lbegin(int level) const {
    assert( level >= 0 );
    // if we dont have this level return empty iterator
    if(level > maxlevel_) return this->template lend<cd,pitype> (level);
    VertexListType & vxList = vertexList_[level];
    if(cd > 1)
    {
      // if vertex list is not up2date, update it
      if(!vxList.up2Date()) vxList.setupVxList(*this,level);
    }
    return ALU3dGridLevelIterator<cd,pitype,const MyType> (*this,vxList,level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU3dGrid<dim, dimworld, elType>::lend(int level) const {
    assert( level >= 0 );
    return ALU3dGridLevelIterator<cd,pitype,const MyType> (*this,level);
  }

  // lbegin methods
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU3dGrid<dim, dimworld, elType>::lbegin(int level) const {
    return this->template lbegin<cd,All_Partition>(level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU3dGrid<dim, dimworld, elType>::lend(int level) const {
    assert( level >= 0 );
    return this->template lend<cd,All_Partition>(level);
  }

  //***********************************************************
  //
  // leaf methods , first all begin methods
  //
  //***********************************************************
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::
  createLeafIteratorBegin(int level) const
  {
    assert( level >= 0 );
    return ALU3dGridLeafIterator<cd, pitype, const MyType> ((*this),level,true);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::leafbegin(int level) const
  {
    return createLeafIteratorBegin<cd, pitype> (level) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin(int level) const {
    return createLeafIteratorBegin<cd, All_Partition> (level) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin() const {
    return createLeafIteratorBegin<codim, pitype> (maxlevel_) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin() const {
    return createLeafIteratorBegin<codim, All_Partition> (maxlevel_) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin(int level) const {
    return createLeafIteratorBegin<0, All_Partition> (level) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin() const {
    return createLeafIteratorBegin<0, All_Partition> (maxlevel_) ;
  }

  //****************************************************************
  //
  // all leaf end methods
  //
  //****************************************************************
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int cd, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::
  createLeafIteratorEnd(int level) const
  {
    assert( level >= 0 );
    return ALU3dGridLeafIterator<cd, pitype, const MyType> ((*this),level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::leafend(int level) const
  {
    return createLeafIteratorEnd <codim, pitype> (level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::leafend(int level) const {
    return createLeafIteratorEnd <codim, All_Partition> (level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::leafend() const {
    return createLeafIteratorEnd <codim, pitype> (maxlevel_);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <int codim>
  inline typename ALU3dGrid<dim, dimworld, elType>::Traits::template Codim<codim>::LeafIterator
  ALU3dGrid<dim, dimworld, elType>::leafend() const {
    return createLeafIteratorEnd <codim, All_Partition> (maxlevel_);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::leafend(int level) const {
    return createLeafIteratorEnd <0, All_Partition> (level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::leafend() const {
    return createLeafIteratorEnd <0,All_Partition> (maxlevel_);
  }

  //*****************************************************************

  // global refine
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim,dimworld, elType>::
  mark(int ref, const typename Traits::template Codim<0>::EntityPointer & ep )
  {
    return this->mark(ref,*ep);
  }

  // global refine
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim,dimworld, elType>::
  mark(int ref, const typename Traits::template Codim<0>::Entity & ep )
  {
    bool marked = (this->getRealImplementation(ep)).mark(ref);
    if(marked)
    {
      if(ref > 0) ++refineMarked_;
      if(ref < 0) ++coarsenMarked_;
    }
    return marked;
  }

  // global refine
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::globalRefine(int numberOfRefines)
  {
    assert( (numberOfRefines + maxLevel()) < MAXL );

    bool ref = false;
    for (int count = numberOfRefines; count>0; count--)
    {
      LeafIteratorType endit  = leafend   ( maxLevel() );
      for(LeafIteratorType it = leafbegin ( maxLevel() ); it != endit; ++it)
      {
        this->mark(1, (*it) );
      }
      ref = this->adapt();
      if(ref) this->postAdapt();
    }

    // important that loadbalance is called on each processor
    // so dont put any if statements arround here
    // this->loadBalance();

    return ref;
  }

  // preprocess grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::preAdapt()
  {
    return (coarsenMarked_ > 0);
  }

  // adapt grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::adapt()
  {
    bool ref = false;
    // if prallel run, then adapt also global id set
    if(globalIdSet_ && comm().size() > 1)
    {
      //std::cout << "Start adapt with globalIdSet prolong \n";
      int defaultChunk = newElementsChunk_;
      int actChunk     = refineEstimate_ * refineMarked_;

      // guess how many new elements we get
      int newElements = std::max( actChunk , defaultChunk );

      globalIdSet_->setChunkSize( newElements );
      //ref = myGrid().adaptWithDataAdaptation(*globalIdSet_); // adapt grid
      ref = myGrid().duneAdapt(*globalIdSet_); // adapt grid
    }
    else
    {
      ref = myGrid().adaptWithoutLoadBalancing();
    }

    if(ref)
    {
      // calcs maxlevel and other extras
      updateStatus();
    }
    return ref;
  }

  // adapt grid
  // --adapt
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <class DofManagerType, class RestrictProlongOperatorType>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  adapt(DofManagerType & dm, RestrictProlongOperatorType & rpo, bool verbose )
  {
    assert( ((verbose) ? (dverb << "ALU3dGrid :: adapt() new method called!\n", 1) : 1 ) );

    typedef typename EntityObject :: ImplementationType EntityImp;
    EntityObject f  ( EntityImp(*this, this->maxLevel()) );
    EntityObject s  ( EntityImp(*this, this->maxLevel()) );

    typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
    typedef CombinedAdaptProlongRestrict < IndexSetRPType,RestrictProlongOperatorType > COType;
    COType tmprpop ( dm.indexSetRPop() , rpo );

    int defaultChunk = newElementsChunk_;
    int actChunk     = refineEstimate_ * refineMarked_;

    // guess how many new elements we get
    int newElements = std::max( actChunk , defaultChunk );

    // reserve memory
    dm.reserveMemory( newElements );

    bool ref = false ;
    if(globalIdSet_ && comm().size() > 1)
    {
      // if global id set exists then include into
      // prolongation process
      ALU3DSPACE AdaptRestrictProlongGlSet<ALU3dGrid<dim, dimworld, elType>,
          COType, GlobalIdSetImp  >
      rp(*this,
         f,this->getRealImplementation(f),
         s,this->getRealImplementation(s),
         tmprpop,
         *globalIdSet_);

      ref = myGrid().duneAdapt(rp); // adapt grid
      if(rp.maxLevel() >= 0) maxlevel_ = rp.maxLevel();
    }
    else
    {
      ALU3DSPACE AdaptRestrictProlongImpl<ALU3dGrid<dim, dimworld, elType>,
          COType >
      rp(*this,
         f,this->getRealImplementation(f),
         s,this->getRealImplementation(s),
         tmprpop);

      ref = myGrid().duneAdapt(rp); // adapt grid
      if(rp.maxLevel() >= 0) maxlevel_ = rp.maxLevel();
    }

    // if new maxlevel was claculated
    assert( ((verbose) ? (dverb << "maxlevel = " << maxlevel_ << "!\n", 1) : 1 ) );

    if(ref)
    {
      updateStatus();
    }

    // check whether we have balance
    dm.dofCompress();

    postAdapt();
    assert( ((verbose) ? (dverb << "ALU3dGrid :: adapt() new method finished!\n", 1) : 1 ) );
    return ref;
  }


  // post process grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::postAdapt()
  {
#if ALU3DGRID_PARALLEL
    {
      // we have to walk over all hierarchcy because during loadBalance
      // we get newly refined elements, which have to be cleared
      int fakeLevel = maxlevel_;
      maxlevel_ = 0;
      for(int l=0; l<= fakeLevel; l++)
      {
        {
          VertexListType & vxList = vertexList_[l];
          typedef ALU3DSPACE ALU3dGridLevelIteratorWrapper<0,Interior_Partition> IteratorType;
          IteratorType w ( *this,vxList,l ,nlinks() ) ;
          for (w.first () ; ! w.done () ; w.next ())
          {
            typedef typename IteratorType :: val_t val_t;
            val_t & item = w.item();
            if(item.first)
            {
              if(item.first->level() > maxlevel_ ) maxlevel_ = item.first->level();
              item.first->resetRefinedTag();
            }
          }
        }
      }

      ALU3DSPACE BSLeafIteratorMaxLevel w ( myGrid() ) ;
      for (w->first () ; ! w->done () ; w->next ())
      {
        if(w->item().level() > maxlevel_ ) maxlevel_ = w->item().level();
        w->item ().resetRefinedTag();

        // note, resetRefinementRequest sets the request to coarsen
        //w->item ().resetRefinementRequest();
      }
    }
#else
    //  if(mpAccess_.nlinks() < 1)
    //#endif
    {
      maxlevel_ = 0;
      ALU3DSPACE BSLeafIteratorMaxLevel w ( myGrid() ) ;
      for (w->first () ; ! w->done () ; w->next ())
      {
        if(w->item().level() > maxlevel_ ) maxlevel_ = w->item().level();
        w->item ().resetRefinedTag();

        // note, resetRefinementRequest sets the request to coarsen
        //w->item ().resetRefinementRequest();
      }
    }
#endif
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::loadBalance()
  {
    if( comm().size() <= 1 ) return false ;
#if ALU3DGRID_PARALLEL
    bool changed = myGrid().duneLoadBalance();
    if(changed)
    {
      // std::cout << "Grid was balanced on p = " << myRank() << std::endl;
      // calculate new maxlevel
      // reset size and things
      myGrid().duneExchangeDynamicState();
      updateStatus();
    }
    return changed;
#else
    return false;
#endif
  }

  // adapt grid
  template <int dim, int dimworld, ALU3dGridElementType elType> template <class DataCollectorType>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  loadBalance(DataCollectorType & dc)
  {
    if( comm().size() <= 1 ) return false ;
#if ALU3DGRID_PARALLEL

    typedef typename EntityObject :: ImplementationType EntityImp;
    EntityObject en     ( EntityImp(*this, this->maxLevel()) );
    EntityObject father ( EntityImp(*this, this->maxLevel()) );
    EntityObject son    ( EntityImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData< ALU3dGrid<dim, dimworld, elType>, DataCollectorType , 0>
    gs(*this,en,this->getRealImplementation(en),dc);

    ALU3DSPACE LoadBalanceRestrictProlongImpl < MyType , DataCollectorType >
    idxop( *this,
           father , this->getRealImplementation(father),
           son    , this->getRealImplementation( son ) ,
           dc );

    // reserve memory for unpacking data
    int defaultChunk = newElementsChunk_;
    int memSize = std::max( idxop.newElements(), defaultChunk );
    dc.reserveMemory ( memSize );

    // call load Balance
    bool changed = myGrid().duneLoadBalance(gs,idxop);

    if(changed)
    {
      dverb << "Grid was balanced on p = " << myRank() << std::endl;
      // calculate new maxlevel
      // reset size and things
    }

    // checken, ob wir das hier wirklich brauchen
    // problem!!!!!!!!!!!!!!!!!!!!!!!
    updateStatus();
    myGrid().duneExchangeData(gs);
    return changed;
#else
    return false;
#endif
  }

  // communicate level data
  template <int dim, int dimworld, ALU3dGridElementType elType> template <class DataHandle>
  inline void ALU3dGrid<dim, dimworld, elType>::
  communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level ) const
  {
    // if only one process, no communication needed
    if( comm().size() <= 1 ) return ;

#if ALU3DGRID_PARALLEL
    const LevelIndexSetImp * levelISet = 0;
    LevelIndexSetImp * newSet = 0;

    if( levelIndexVec_[level] == 0 )
    {
      newSet = new LevelIndexSetImp ( *this , level );
      levelISet = newSet;
    }
    else
      levelISet = levelIndexVec_[level] ;

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim>::Entity> VertexObject;
    typedef typename Traits :: template Codim<dim>::Entity::ImplementationType VertexImp;
    VertexObject vx( VertexImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandle, dim>
    vertexData(*this,vx,this->getRealImplementation(vx),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim-1>::Entity> EdgeObject;
    typedef typename Traits::template Codim<dim-1>::Entity::ImplementationType EdgeImp;
    EdgeObject edge( EdgeImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandle, dim-1>
    edgeData(*this,edge,this->getRealImplementation(edge),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<1>::Entity> FaceObject;
    typedef typename Traits::template Codim<1>::Entity::ImplementationType FaceImp;
    FaceObject face( FaceImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandle, 1>
    faceData(*this,face,this->getRealImplementation(face),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<0>::Entity> ElementObject;
    typedef typename Traits::template Codim<0>::Entity::ImplementationType ElementImp;
    ElementObject element( ElementImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandle, 0>
    elementData(*this,element,this->getRealImplementation(element),data,*levelISet,level);

    doCommunication( vertexData, edgeData, faceData, elementData , iftype, dir );

    if(newSet) delete newSet;
#endif
  }

  // communicate data
  template <int dim, int dimworld, ALU3dGridElementType elType> template <class DataHandle>
  inline void ALU3dGrid<dim, dimworld, elType>::
  communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
  {
    // if only one process, no communication needed
    if( comm().size() <= 1 ) return ;

#if ALU3DGRID_PARALLEL

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim>::Entity> VertexObject;
    typedef typename Traits :: template Codim<dim>::Entity::ImplementationType VertexImp;
    VertexObject vx( VertexImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandle, dim>
    vertexData(*this,vx,this->getRealImplementation(vx),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim-1>::Entity> EdgeObject;
    typedef typename Traits::template Codim<dim-1>::Entity::ImplementationType EdgeImp;
    EdgeObject edge( EdgeImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandle, dim-1>
    edgeData(*this,edge,this->getRealImplementation(edge),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<1>::Entity> FaceObject;
    typedef typename Traits::template Codim<1>::Entity::ImplementationType FaceImp;
    FaceObject face( FaceImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandle, 1>
    faceData(*this,face,this->getRealImplementation(face),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<0>::Entity> ElementObject;
    typedef typename Traits::template Codim<0>::Entity::ImplementationType ElementImp;
    ElementObject element( ElementImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandle, 0>
    elementData(*this,element,this->getRealImplementation(element),data);

    doCommunication( vertexData, edgeData, faceData, elementData , iftype, dir );
#endif
  }

  // communicate data
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::
  doCommunication ( GatherScatterType & vertexData,
                    GatherScatterType & edgeData,
                    GatherScatterType & faceData,
                    GatherScatterType & elementData,
                    InterfaceType iftype, CommunicationDirection dir) const
  {
    // check interface types
    if( (iftype == Overlap_OverlapFront_Interface) ||
        (iftype == Overlap_All_Interface) )
    {
      dverb << "ALUGrid contains no overlap, therefore no communication for\n";
      dverb << "Overlap_OverlapFront_Interface or Overlap_All_Interface interfaces! \n";
      return ;
    }
#if ALU3DGRID_PARALLEL
    // communication from border to border
    if( iftype == InteriorBorder_InteriorBorder_Interface )
    {
      myGrid().borderBorderCommunication(vertexData,edgeData,faceData,elementData);
      return ;
    }

    // communication from interior to ghost including border
    if( iftype == InteriorBorder_All_Interface )
    {
      if( dir == ForwardCommunication )
      {
        myGrid().interiorGhostCommunication(vertexData,edgeData,faceData,elementData);
        return ;
      }
      // reverse communiction interface (here All_InteriorBorder)
      if( dir == BackwardCommunication )
      {
        myGrid().ghostInteriorCommunication(vertexData,edgeData,faceData,elementData);
        return ;
      }
    }

    // communication from interior to ghost including border
    if( iftype == All_All_Interface )
    {
      myGrid().allAllCommunication(vertexData,edgeData,faceData,elementData);
      return ;
    }

    DUNE_THROW(GridError,"wrong set of parameters in ALUGrid::doCommunication");
#endif
  }


  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  writeGrid_Ascii(const std::string filename, alu3d_ctype time ) const
  {
    ALU3DSPACE GitterImplType & mygrd = myGrid();
    std::fstream file ( filename.c_str() , std::ios::out);
    if(file)
    {
      typedef typename ALU3dImplTraits<elType> :: BNDFaceType BNDFaceType;
      typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;
      typedef typename ALU3dImplTraits<elType> :: HasFaceType HasFaceType;

      file << "!" << elType2Name( elType ) << std::endl;
      {
        ALU3DSPACE LeafIterator < ALU3DSPACE VertexType > vx (mygrd) ;
        file << std::endl;

        // write coordinates of the vertices
        int vxsize = vx->size();
        file << vxsize << std::endl;
        Array < double[3] > vxvec ( vxsize );

        for( vx->first(); !vx->done() ; vx->next() )
        {
          const double (&p)[3] = vx->item().Point();
          int vxidx = vx->item().getIndex();
          double (&v)[3] = vxvec[vxidx];
          for(int i=0; i<3; i++) v[i] = p[i];
        }

        for(int i=0; i<vxsize; i++)
        {
          file << vxvec[i][0] << " " << vxvec[i][1] << " " << vxvec[i][2] << std::endl;
        }
      }

      file << std::endl;
      // write element vertices
      {
        const int novx = (elType == tetra) ? 4 : 8;
        ALU3DSPACE LeafIterator < ALU3DSPACE HElementType > el (mygrd) ;
        file << el->size() << std::endl;
        for( el->first(); !el->done() ; el->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (el->item());
          for(int i=0; i<novx; i++)
          {
            const int vxnum = item.myvertex(i)->getIndex();
            file << vxnum << " ";
          }
          file << std::endl;
        }
      }

      // write boundary faces
      {
        file << std::endl;
        const int nofaces  = (elType == tetra) ? 4 : 6;
        int bndfaces = 0;
        ALU3DSPACE LeafIterator < ALU3DSPACE HElementType > el (mygrd) ;
        for( el->first(); !el->done() ; el->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (el->item());
          for(int i=0; i<nofaces; i++)
          {
            std::pair < HasFaceType * , int > nbpair = item.myneighbour(i);
            if(nbpair.first->isboundary())
            {
              bndfaces++;
            }
          }
        }
        file << bndfaces << std::endl;
      }
      // write boundary faces
      {
        const int bndvxnum = (elType == tetra) ? 3 : 4;
        const int nofaces  = (elType == tetra) ? 4 : 6;
        ALU3DSPACE LeafIterator < ALU3DSPACE HElementType > el (mygrd) ;
        for( el->first(); !el->done() ; el->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (el->item());
          for(int i=0; i<nofaces; i++)
          {
            std::pair < HasFaceType * , int > nbpair = item.myneighbour(i);
            if(nbpair.first->isboundary())
            {
              BNDFaceType * face = static_cast<BNDFaceType *> (nbpair.first);
              file << -(face->bndtype()) << " " << bndvxnum << " ";
              for(int j=0; j<bndvxnum; j++)
              {
                int vxnum = face->myvertex(0,j)->getIndex();
                file << vxnum << " ";
              }
              file << std::endl;
            }
          }
        }
      }

      {
        ALU3DSPACE LeafIterator < ALU3DSPACE VertexType > vx (mygrd) ;
        file << std::endl;

        // write coordinates of the vertices
        int vxnum = 0;
        for( vx->first(); !vx->done() ; vx->next() )
        {
          file << vxnum << " -1" << std::endl;
          vxnum++;
        }
      }
    }
    return true;
  }


  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <GrapeIOFileFormatType ftype>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  writeGrid(const std::string filename, alu3d_ctype time ) const
  {
    switch(ftype)
    {
    case xdr  : return writeGrid_Xdr(filename,time);
    case ascii : return writeGrid_Ascii(filename,time);
    default : derr << "Wrong file type in writeGrid method~ \n";
    }
    return false;
  }


  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  writeGrid_Xdr(const std::string filename, alu3d_ctype time ) const
  {
    ALU3DSPACE GitterImplType & mygrd = myGrid();
    mygrd.duneBackup(filename.c_str());

    // write time and maxlevel
    {
      std::string extraName(filename);
      extraName += ".extra";
      std::ofstream out (extraName.c_str());
      if(out)
      {
        out.precision (16);
        out << time << " ";
        out << maxlevel_ << " ";
        out.close();
      }
      else
      {
        derr << "ALU3dGrid::writeGrid: couldn't open <" << extraName << ">! \n";
      }
    }
    return true;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <GrapeIOFileFormatType ftype>
  inline bool ALU3dGrid<dim,dimworld, elType>::
  readGrid( const std::string filename, alu3d_ctype & time )
  {
    {
      typedef std::ostringstream StreamType;
      std::string mName(filename);
      mName += ".macro";
      const char * macroName = mName.c_str();

      { //check if file exists
        std::ifstream check ( macroName );
        if( !check )
          DUNE_THROW(GridError,"cannot read file " << macroName << "\n");
        check.close();
      }

      mygrid_ = new ALU3DSPACE GitterImplType (macroName
#if ALU3DGRID_PARALLEL
                                               , mpAccess_
#endif
                                               );
    }


    assert(mygrid_ != 0);

    // check for element type
    this->checkMacroGrid ();

    myGrid().duneRestore(filename.c_str());

    {
      std::string extraName (filename);
      extraName += ".extra";
      std::ifstream in (extraName.c_str());
      if(in)
      {
        in.precision (16);
        in  >> time;
        in  >> maxlevel_;
        in.close();
      }
      else
      {
        derr << "ALU3dGrid::readGrid: couldn't open <" << extraName << ">! \n";
      }
    }

    // calculate new maxlevel
    // calculate indices
    updateStatus();

    return true;
  }

  // return Grid type
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline GridIdentifier ALU3dGrid<dim, dimworld, elType>::type () const
  {
    return ALU3dGrid_Id;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::checkMacroGrid()
  {
    typedef ALU3DSPACE LevelIterator < ALU3DSPACE HElementType > IteratorType;
    IteratorType w( this->myGrid(), 0 );
    for (w->first () ; ! w->done () ; w->next ())
    {
      ALU3dGridElementType type = (ALU3dGridElementType) w->item().type();
      if( type != elType )
      {
        derr << "\nERROR: " << elType2Name(elType) << " Grid tries to read a ";
        derr << elType2Name(type) << " macro grid file! \n\n";
        assert(type == elType);
        abort();
      }
    }
  }

  inline const char * elType2Name( ALU3dGridElementType elType )
  {
    switch( elType )
    {
    case tetra  : return "Tetraeder";
    case hexa   : return "Hexaeder";
    case mixed  : return "Mixed";
    default     : return "Error";
    }
  }

} // end namespace Dune
