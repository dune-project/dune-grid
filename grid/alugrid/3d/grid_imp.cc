// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Dune includes
#include <dune/common/stdstreams.hh>
#include <dune/common/interfaces.hh>

// Local includes
#include "entity.hh"
#include "iterator.hh"
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
    std::vector<int> visited_(vxsize);

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

    ElementLevelIteratorType it ( grid, level, grid.nlinks() );

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

  template <class GridType >
  inline void ALU3dGridLeafVertexList ::
  setupVxList(const GridType & grid)
  {
    // iterates over grid elements of given level and adds all vertices to
    // given list
    enum { codim = 3 };

    VertexListType & vxList = vertexList_;
    size_t vxsize = grid.hierarchicIndexSet().size(codim);
    if( vxList.capacity() < vxsize) vxList.reserve(vxsize);
    vxList.resize(vxsize);

    for(size_t i=0; i<vxsize; ++i)
    {
      ItemType & vx = vxList[i];
      vx.first  = 0;
      vx.second = -1;
    }

    const ALU3dGridElementType elType = GridType:: elementType;

    typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper<0,Dune::All_Partition>  ElementIteratorType;
    typedef typename ElementIteratorType :: val_t val_t;

    typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;
    typedef ALU3DSPACE VertexType VertexType;

    enum { nVx = ElementTopologyMapping < elType > :: numVertices };

    ElementIteratorType it ( grid, grid.maxLevel() , grid.nlinks() );

#ifndef NDEBUG
    int count = 0;
#endif
    for( it.first(); !it.done() ; it.next())
    {
      val_t & item = it.item();

      IMPLElementType * elem = 0;
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first);
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first);

      assert( elem );
      int level = elem->level();

      for(int i=0; i<nVx; ++i)
      {
        VertexType * vx = elem->myvertex(i);
        assert( vx );

        // insert only interior and border vertices
        if( vx->isGhost() ) continue;

        const int idx = vx->getIndex();
        ItemType & vxpair = vxList[idx];
        if( vxpair.first == 0 )
        {
          vxpair.first  = vx;
          vxpair.second = level;
#ifndef NDEBUG
          ++ count;
#endif
        }
        // always store max level of vertex as grdi definition says
        else
        {
          // set the max level for each vertex, see Grid definition
          if (vxpair.second < level) vxpair.second = level;
        }
      }
    }

    //std::cout << count << "c | s " << vxList.size() << "\n";
    // make sure that the found number of vertices equals to stored ones
    //assert( count == (int)vxList.size() );
    up2Date_ = true;
  }


  //--Grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::
  ALU3dGrid(const std::string macroTriangFilename,
#if ALU3DGRID_PARALLEL
            const MPI_Comm mpiComm,
#endif
            const DuneBoundaryProjectionType* bndPrj
            )
    : mygrid_ (0)
#if ALU3DGRID_PARALLEL
      , mpAccess_(mpiComm)
      , ccobj_(mpiComm)
#endif
      , maxlevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_ (*this)
      , globalIdSet_(0), localIdSet_(*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , sizeCache_ (0)
      , lockPostAdapt_(false)
      , vertexProjection_( (bndPrj) ? new ALUGridBoundaryProjectionType( *bndPrj ) : 0 )
  {
    makeGeomTypes();

    // check macro grid file for keyword
    this->checkMacroGridFile (macroTriangFilename);

#if ALU3DGRID_PARALLEL == 0
    // in serial version call empty constructor when no file given
    if( macroTriangFilename == "" )
    {
      mygrid_ = new ALU3DSPACE GitterImplType ();
    }
    else
#endif
    {
      mygrid_ = new ALU3DSPACE GitterImplType (macroTriangFilename.c_str()
#if ALU3DGRID_PARALLEL
                                               , mpAccess_ ,
#endif
#ifdef ALUGRID_VERTEX_PROJECTION
                                               , vertexProjection_ // only for newer versions
#endif
                                               );
    }

    assert(mygrid_ != 0);

#if ALU3DGRID_PARALLEL
    dverb << "************************************************\n";
    dverb << "Created grid on p=" << mpAccess_.myrank() << "\n";
    dverb << "************************************************\n";
#endif
    this->checkMacroGrid ();

    postAdapt();
    calcExtras();
  } // end Constructor

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType> & ALU3dGrid<dim, dimworld, elType>::operator = (const ALU3dGrid<dim, dimworld, elType> & g)
  {
    DUNE_THROW(GridError,"Do not use assignment operator of ALU3dGrid! \n");
    return (*this);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::~ALU3dGrid()
  {
    delete vertexProjection_;

    for(unsigned int i=0; i<levelIndexVec_.size(); i++) delete levelIndexVec_[i];
    delete globalIdSet_; globalIdSet_ = 0;
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;

    /*
       if(myGrid().container().iterators_attached())
       {
       dwarn << "WRANING: There still exists instances of iterators giving access to this grid which is about to be removed! in: " << __FILE__ << " line: " << __LINE__ << std::endl;
       }
     */
    delete mygrid_; mygrid_ = 0;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::size(int level, int codim) const
  {
    // if we dont have this level return 0
    if( (level > maxlevel_) || (level < 0) ) return 0;

    assert( codim >= 0);
    assert( codim < dim+1 );

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

  // --size
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

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::ghostSize(int codim) const
  {
#if ALU3DGRID_PARALLEL
    assert( codim >= 0 );
    assert( codim < dim +1 );

    assert( sizeCache_ );
    return sizeCache_->ghostSize(codim);
#else
    return 0;
#endif
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::ghostSize(int level, int codim) const
  {
#if ALU3DGRID_PARALLEL
    assert( codim >= 0 );
    assert( codim < dim +1 );
    assert( level >= 0);

    assert( sizeCache_ );
    return sizeCache_->ghostSize(level,codim);
#else
    return 0;
#endif
  }

  // calc all necessary things that might have changed
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::setMaxLevel(int mxl)
  {
    maxlevel_ = std::max(maxlevel_,mxl);
  }

  // calc all necessary things that might have changed
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::updateStatus()
  {
    calcMaxLevel();
    calcExtras();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::calcMaxLevel()
  {
    // old fashioned way
    int testMaxLevel = 0;
    typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper<0,All_Partition> IteratorType;
    IteratorType w (*this, maxLevel(), nlinks() );

    typedef typename IteratorType :: val_t val_t ;
    typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;

    for (w.first () ; ! w.done () ; w.next ())
    {
      val_t & item = w.item();

      IMPLElementType * elem = 0;
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first);
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first);

      assert( elem );

      int level = elem->level();
      if(level > testMaxLevel) testMaxLevel = level;
    }
    maxlevel_ = testMaxLevel;
  }

  // --calcExtras
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::calcExtras()
  {
    if(sizeCache_) delete sizeCache_;
    bool isSimplex = (elType == tetra) ? true : false;
    sizeCache_ = new SizeCacheType (*this,isSimplex,!isSimplex,true);

    // unset up2date before recalculating the index sets,
    // becasue they will use this feature
    leafVertexList_.unsetUp2Date();
    for(size_t i=0; i<MAXL; ++i)
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

    // update all index set that are already in use
    for(size_t i=0; i<levelIndexVec_.size(); ++i)
    {
      if(levelIndexVec_[i]) (*(levelIndexVec_[i])).calcNewIndex();
    }
    if(leafIndexSet_) leafIndexSet_->calcNewIndex();

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

    return ALU3dGridLevelIterator<cd,pitype,const MyType> (*this,level,true);
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

  // mark given entity
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

  // get Mark of given entity
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim,dimworld, elType>::
  getMark(const typename Traits::template Codim<0>::Entity & en) const
  {
    return this->getRealImplementation(en).getMark();
  }

  // global refine
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::globalRefine ( int refCount )
  {
    assert( (refCount + maxLevel()) < MAXL );

    for( int count = refCount; count > 0; --count )
    {
      const LeafIteratorType end = leafend();
      for( LeafIteratorType it = leafbegin(); it != end; ++it )
        mark( 1, *it );
      const bool refined = adapt();
      if( refined )
        postAdapt();
    }
  }


  // global refine
  template< int dim, int dimworld, ALU3dGridElementType elType >
  template< class GridImp, class DataHandle >
  inline void ALU3dGrid< dim, dimworld, elType >
  ::globalRefine ( int refCount, AdaptDataHandleInterface< GridImp, DataHandle > &handle )
  {
    assert( (refCount + maxLevel()) < MAXL );

    for( int count = refCount; count > 0; --count )
    {
      const LeafIteratorType end = leafend();
      for( LeafIteratorType it = leafbegin(); it != end; ++it )
        mark( 1 , *it );
      adapt( handle );
    }
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

    if( lockPostAdapt_ == true )
    {
      DUNE_THROW(InvalidStateException,"Make sure that postAdapt is called after adapt was called and returned true!");
    }

    bool mightCoarse = preAdapt();
    // if prallel run, then adapt also global id set
    if(globalIdSet_)
    {
      //std::cout << "Start adapt with globalIdSet prolong \n";
      int defaultChunk = newElementsChunk_;
      int actChunk     = refineEstimate_ * refineMarked_;

      // guess how many new elements we get
      int newElements = std::max( actChunk , defaultChunk );

      globalIdSet_->setChunkSize( newElements );
      ref = myGrid().duneAdapt(*globalIdSet_); // adapt grid
    }
    else
    {
      ref = myGrid().adaptWithoutLoadBalancing();
    }

    // in parallel this is different
    if( this->comm().size() == 1 )
    {
      ref = ref && refineMarked_ > 0;
    }

    if(ref || mightCoarse)
    {
      // calcs maxlevel and other extras
      updateStatus();

      // notify that postAdapt must be called
      lockPostAdapt_ = true;
    }
    return ref;
  }


  // adapt grid
  // --adapt
  template< int dim, int dimworld, ALU3dGridElementType elType >
  template< class GridImp, class DataHandle >
  inline bool ALU3dGrid< dim, dimworld, elType >
  ::adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle )
  {
    typedef AdaptDataHandleInterface< GridImp, DataHandle > AdaptDataHandle;

    typedef typename EntityObject::ImplementationType EntityImp;
    EntityObject father( EntityImp( *this, this->maxLevel() ) );
    EntityObject son( EntityImp( *this, this->maxLevel() ) );

    int defaultChunk = newElementsChunk_;
    int actChunk     = refineEstimate_ * refineMarked_;

    // guess how many new elements we get
    int newElements = std::max( actChunk , defaultChunk );

    // true if at least one element was marked for coarsening
    bool mightCoarse = preAdapt();
    // reserve memory
    handle.preAdapt( newElements );

    bool refined = false ;
    if(globalIdSet_)
    {
      // if global id set exists then include into
      // prolongation process
      ALU3DSPACE AdaptRestrictProlongGlSet< MyType, AdaptDataHandle, GlobalIdSetImp >
      rp(*this,
         father,this->getRealImplementation(father),
         son,   this->getRealImplementation(son),
         handle,
         *globalIdSet_);

      refined = myGrid().duneAdapt(rp); // adapt grid
    }
    else
    {
      ALU3DSPACE AdaptRestrictProlongImpl< MyType, AdaptDataHandle >
      rp(*this,
         father,this->getRealImplementation(father),
         son,   this->getRealImplementation(son),
         handle);

      refined = myGrid().duneAdapt(rp); // adapt grid
    }

    if(refined || mightCoarse)
    {
      // only calc extras and skip maxLevel calculation, because of
      // refinement maxLevel was calculated already
      updateStatus();

      // no need to call postAdapt here, because markers
      // are cleand during refinement callback
    }

    // check whether we have balance
    handle.postAdapt();

    // here postAdapt is not called, because
    // reset of refinedTag is done in preCoarsening and postRefinement
    // methods of datahandle (see datahandle.hh)

    return refined;
  }

  // post process grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::postAdapt()
  {
    {
      // old fashioned way
      typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper<0,All_Partition> IteratorType;
      IteratorType w (*this, maxLevel(), nlinks() );

      typedef typename IteratorType :: val_t val_t ;
      typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;

      for (w.first () ; ! w.done () ; w.next ())
      {
        val_t & item = w.item();

        IMPLElementType * elem = 0;
        if( item.first )
          elem = static_cast<IMPLElementType *> (item.first);
        else if( item.second )
          elem = static_cast<IMPLElementType *> (item.second->getGhost().first);

        assert( elem );
        elem->resetRefinedTag();
      }
    }

    // make that postAdapt has been called
    lockPostAdapt_ = false;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline bool ALU3dGrid<dim, dimworld, elType>::loadBalance()
  {
    if( comm().size() <= 1 ) return false ;
#if ALU3DGRID_PARALLEL
    bool changed = myGrid().duneLoadBalance();
    if(changed)
    {
      // some exchanges on ALUGrid side
      myGrid().duneExchangeDynamicState();

      // calculate new maxlevel
      // reset size and things
      updateStatus();

      // build new Id Set. Only do that after updateStatus, because here
      // the item lists are needed
      if(globalIdSet_) globalIdSet_->updateIdSet();

      // unset all leaf markers
      postAdapt();
    }

    return changed;
#else
    return false;
#endif
  }

  // load balance grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <class DataHandle>
  inline bool ALU3dGrid<dim, dimworld, elType>::
  loadBalance(DataHandle & data)
  {
    if( comm().size() <= 1 ) return false ;
#if ALU3DGRID_PARALLEL
    typedef typename EntityObject :: ImplementationType EntityImp;
    EntityObject en     ( EntityImp(*this, this->maxLevel()) );
    EntityObject father ( EntityImp(*this, this->maxLevel()) );
    EntityObject son    ( EntityImp(*this, this->maxLevel()) );

    typedef ALU3DSPACE LoadBalanceElementCount<ThisType,DataHandle> LDBElCountType;

    // elCount is the adaption restPro operator used during the refinement
    // cause be creating new elements on processors
    LDBElCountType elCount(*this,
                           father,this->getRealImplementation(father),
                           son,this->getRealImplementation(son),
                           data);

    ALU3DSPACE GatherScatterLoadBalance< ALU3dGrid<dim, dimworld, elType>,
        DataHandle, LDBElCountType >
    gs(*this,en,this->getRealImplementation(en),data,elCount);

    // call load Balance
    bool changed = myGrid().duneLoadBalance(gs,elCount);

    if(changed)
    {
      // exchange some data for internal useage
      myGrid().duneExchangeDynamicState();

      // calculate new maxlevel
      // reset size and things
      updateStatus();

      // build new Id Set. Only do that after updateStatus, because here
      // the item lists are needed
      if(globalIdSet_) globalIdSet_->updateIdSet();

      // compress data , wrapper for dof manager
      gs.compress();

      postAdapt();
    }
    return changed;
#else
    return false;
#endif
  }

  // communicate level data
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <class DataHandleImp,class DataType>
  inline void ALU3dGrid<dim, dimworld, elType>::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir, int level ) const
  {
    // if only one process, no communication needed
    if( comm().size() <= 1 ) return ;

#if ALU3DGRID_PARALLEL
    typedef CommDataHandleIF<DataHandleImp,DataType> DataHandleType;

    // for level communication the level index set is needed.
    // fi non-existent, then create for communicaton
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
    typedef typename VertexObject::ImplementationType VertexImp;
    VertexObject vx( VertexImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandleType, dim>
    vertexData(*this,vx,this->getRealImplementation(vx),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim-1>::Entity> EdgeObject;
    typedef typename EdgeObject::ImplementationType EdgeImp;
    EdgeObject edge( EdgeImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandleType, dim-1>
    edgeData(*this,edge,this->getRealImplementation(edge),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<1>::Entity> FaceObject;
    typedef typename FaceObject::ImplementationType FaceImp;
    FaceObject face( FaceImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandleType, 1>
    faceData(*this,face,this->getRealImplementation(face),data,*levelISet,level);

    typedef MakeableInterfaceObject<typename Traits::template Codim<0>::Entity> ElementObject;
    typedef typename ElementObject::ImplementationType ElementImp;
    ElementObject element( ElementImp(*this, level ) );

    ALU3DSPACE GatherScatterLevelData < ThisType, DataHandleType, 0>
    elementData(*this,element,this->getRealImplementation(element),data,*levelISet,level);

    doCommunication( vertexData, edgeData, faceData, elementData , iftype, dir );

    if(newSet) delete newSet;
#endif
  }

  // communicate data
  template <int dim, int dimworld, ALU3dGridElementType elType>
  template <class DataHandleImp, class DataType>
  inline void ALU3dGrid<dim, dimworld, elType>::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir) const
  {
    // if only one process, no communication needed
    if( comm().size() <= 1 ) return ;

#if ALU3DGRID_PARALLEL
    typedef CommDataHandleIF<DataHandleImp,DataType> DataHandleType;

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim>::Entity> VertexObject;
    typedef typename VertexObject::ImplementationType VertexImp;
    VertexObject vx( VertexImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandleType, dim>
    vertexData(*this,vx,this->getRealImplementation(vx),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<dim-1>::Entity> EdgeObject;
    typedef typename EdgeObject::ImplementationType EdgeImp;
    EdgeObject edge( EdgeImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandleType, dim-1>
    edgeData(*this,edge,this->getRealImplementation(edge),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<1>::Entity> FaceObject;
    typedef typename FaceObject::ImplementationType FaceImp;
    FaceObject face( FaceImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandleType, 1>
    faceData(*this,face,this->getRealImplementation(face),data);

    typedef MakeableInterfaceObject<typename Traits::template Codim<0>::Entity> ElementObject;
    typedef typename ElementObject::ImplementationType ElementImp;
    ElementObject element( ElementImp(*this, this->maxLevel()) );

    ALU3DSPACE GatherScatterLeafData < ThisType, DataHandleType, 0>
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
  writeGrid_Ascii(const std::string filename, alu3d_ctype time , bool scientific  ) const
  {
#if ALU3DGRID_PARALLEL
    // write ascii only works for serial grids at the moment
    if ( this->comm().size() > 1 )
    {
      DUNE_THROW(GridError,"ALU3dGrid::writeGrid_Ascii not implemented for parallel grids!");
    }
#endif

    ALU3DSPACE GitterImplType & mygrd = myGrid();
    std::ofstream file ( filename.c_str() );
    if(file)
    {
      typedef typename ALU3dImplTraits<elType> :: BNDFaceType BNDFaceType;
      typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;
      typedef typename ALU3dImplTraits<elType> :: HasFaceType HasFaceType;
      typedef typename ALU3dImplTraits<elType> :: GEOVertexType GEOVertexType;

      ALU3DSPACE LeafIterator < ALU3DSPACE HElementType > leafElements (mygrd) ;
      file << "!" << elType2Name( elType ) << " Elements = " << leafElements->size() << std::endl;

      ALU3DSPACE LeafIterator < ALU3DSPACE VertexType > leafVertices (mygrd) ;
      {
        file << std::endl;

        if( scientific )
        {
          // use scientific mode
          file << std::scientific;
        }

        // write coordinates of the vertices
        int vxsize = leafVertices->size();
        file << vxsize << std::endl;

        typedef double ShortVecType[3];
        ShortVecType * vxvec = new ShortVecType [vxsize];
        assert( vxvec );

        for( leafVertices->first(); !leafVertices->done() ; leafVertices->next() )
        {
          const GEOVertexType & vertex =
            static_cast<GEOVertexType &> (leafVertices->item());
          const double (&p)[3] = vertex.Point();
          int vxidx = vertex.getIndex();
          double (&v)[3] = vxvec[vxidx];
          for(int i=0; i<3; ++i) v[i] = p[i];
        }

        for(int i=0; i<vxsize; ++i)
        {
          file << vxvec[i][0] << " " << vxvec[i][1] << " " << vxvec[i][2] << std::endl;
        }

        delete [] vxvec;
      }

      file << std::endl;
      // write element vertices
      {
        const int novx = (elType == tetra) ? 4 : 8;
        file << leafElements->size() << std::endl;
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<novx; ++i)
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
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<nofaces; ++i)
          {
            std::pair < HasFaceType * , int > nbpair = item.myneighbour(i);
            if(nbpair.first->isboundary())
            {
              ++bndfaces;
            }
          }
        }
        file << bndfaces << std::endl;
      }
      // write boundary faces
      {
        const int bndvxnum = (elType == tetra) ? 3 : 4;
        const int nofaces  = (elType == tetra) ? 4 : 6;
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<nofaces; ++i)
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
        file << std::endl;

        // write coordinates of the vertices
        int vxnum = 0;
        for( leafVertices->first(); !leafVertices->done() ; leafVertices->next() )
        {
          file << vxnum << " -1" << std::endl;
          ++vxnum;
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
        out << std::scientific << time << " ";
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
        if( !check
            // only abort on rank 0
            // on all other ranks this can be empty
            && comm().rank() == 0
            )
          DUNE_THROW(GridError,"cannot read file " << macroName << "\n");
        check.close();
      }

      // if grid exists delete first
      if( mygrid_ ) delete mygrid_;
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
        in  >> std::scientific >> time;
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

    // reset refinement markers
    postAdapt();

    // send time from proc 0 to all in case that some grids are empty
    comm().broadcast(&time,1,0);

    return true;
  }

  // return Grid type
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline std::string ALU3dGrid<dim, dimworld, elType>::name () const
  {
    if(elType == hexa)
      return "ALUCubeGrid";
    assert( elType == tetra );
    return "ALUSimplexGrid";
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline void ALU3dGrid<dim, dimworld, elType>::checkMacroGridFile(const std::string filename)
  {
    if(filename == "") return;

    std::ifstream file(filename.c_str());
    if(!file)
    {
      std::cerr << "Couldn't open file '" << filename <<"' !" << std::endl;
      DUNE_THROW(IOError,"Couldn't open file '" << filename <<"' !");
    }

    const std::string aluid((elType == tetra) ? "!Tetrahedra" : "!Hexahedra");
    const std::string oldAluId((elType == tetra) ? "!Tetraeder" : "!Hexaeder");
    std::string idline;
    std::getline(file,idline);
    std::stringstream idstream(idline);
    std::string id;
    idstream >> id;

    if(id == aluid )
    {
      return;
    }
    else if ( id == oldAluId )
    {
      derr << "\nKeyword '" << oldAluId << "' is deprecated! Change it to '" << aluid << "' in file '" << filename<<"'! \n";
      return ;
    }
    else
    {
      std::cerr << "Delivered file '"<<filename<<"' does not contain keyword '"
                << aluid << "'. Found id '" <<id<< "'. Check the macro grid file! Bye." << std::endl;
      DUNE_THROW(IOError,"Wrong file format! ");
    }
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
        DUNE_THROW(GridError,"\nERROR: " << elType2Name(elType) << " Grid tries to read a " << elType2Name(type) << " macro grid file! ");
      }
    }
  }

  inline const char * elType2Name( ALU3dGridElementType elType )
  {
    switch( elType )
    {
    case tetra  : return "Tetrahedra";
    case hexa   : return "Hexahedra";
    case mixed  : return "Mixed";
    default     : return "Error";
    }
  }

} // end namespace Dune
