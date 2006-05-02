// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALU2DGRID_IMP_CC
#define DUNE_ALU2DGRID_IMP_CC

namespace Dune {

  /*
     template <class GridType >
     inline void ALU2dGridVertexList ::
     setupVxList(const GridType & grid, int level)
     {
      // iterates over grid elements of given level and adds all vertices to
      // given list
      enum { codim = 3 };

      VertexListType & vxList = vertexList_;

      unsigned int vxsize = grid.hierarchicIndexSet().size(codim);
      if( vxList.size() < vxsize ) vxList.resize(vxsize);
      for(unsigned int i=0; i<vxsize; i++) vxList[i] = 0;

      const ALU3dGridElementType elType = GridType:: elementType;
      typedef ALU3DSPACE LevelIterator < ALU3DSPACE HElementType > IteratorType;
      typedef typename ALU3dImplTraits<elType> :: IMPLElementType IMPLElementType;
      typedef ALU3DSPACE VertexType VertexType;
      enum { nVx = ElementTopologyMapping < elType > :: numVertices };

      IteratorType it (const_cast<GridType &> (grid).myGrid(),level);
      for( it->first(); !it->done() ; it->next())
      {
        IMPLElementType & elem = static_cast<IMPLElementType &> (it->item());
        for(int i=0; i<nVx; i++)
        {
          VertexType * vx = elem.myvertex(i);
          vxList[vx->getIndex()] = vx;
        }
      }
     }
   */

  //--Grid

  //! Iterator to first entity of given codim on level
  template <int dim, int dimworld>
  template<int cd, PartitionIteratorType pitype>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU2dGrid<dim, dimworld> :: lbegin (int level) const {
    return ALU2dGridLevelIterator<cd, pitype, const ThisType>(*this, level, false);
  }

  //! one past the end on this level
  template <int dim, int dimworld>
  template<int cd, PartitionIteratorType pitype>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU2dGrid<dim, dimworld> :: lend (int level) const {
    return ALU2dGridLevelIterator<cd, pitype, const ThisType>(*this, level, true);
  }

  //! Iterator to first entity of given codim on level
  template <int dim, int dimworld>
  template<int cd>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU2dGrid<dim, dimworld>::lbegin (int level) const {
    return ALU2dGridLevelIterator<cd, All_Partition, const ThisType>(*this, level, false);
  }

  //! one past the end on this level
  template <int dim, int dimworld>
  template<int cd>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU2dGrid<dim, dimworld>::lend (int level) const {
    return ALU2dGridLevelIterator<cd, All_Partition, const ThisType>(*this, level, true);
  }

  //! Iterator to first entity of codim 0 on level
  template <int dim, int dimworld>
  inline typename ALU2dGrid<dim, dimworld>::LevelIteratorType ALU2dGrid<dim, dimworld>::lbegin (int level) const {
    return LevelIteratorImp(*this, level, false);
  }

  //! last entity of codim 0 on level
  template <int dim, int dimworld>
  inline typename ALU2dGrid<dim, dimworld>::LevelIteratorType ALU2dGrid<dim, dimworld> :: lend (int level) const {
    return LevelIteratorImp(*this, level, true);
  }

  //! General definiton for a leaf iterator
  template <int dim, int dimworld>
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU2dGrid<dim, dimworld> :: leafbegin() const {
    return ALU2dGridLeafIterator<codim, pitype, const ThisType> (*this, false);
  }

  //! General definition for an end iterator on leaf level
  template <int dim, int dimworld>
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU2dGrid<dim, dimworld>:: leafend() const {
    return ALU2dGridLeafIterator<codim, pitype, const ThisType> (*this, true);
  }

  //! General definiton for a leaf iterator
  template <int dim, int dimworld>
  template <int codim>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<codim>::LeafIterator
  ALU2dGrid<dim, dimworld>:: leafbegin() const {
    return ALU2dGridLeafIterator<codim, All_Partition, const ThisType> (*this, false);
  }

  //! General definition for an end iterator on leaf level
  template <int dim, int dimworld>
  template <int codim>
  inline typename ALU2dGrid<dim, dimworld>::Traits::template Codim<codim>::LeafIterator
  ALU2dGrid<dim, dimworld>:: leafend() const {
    return ALU2dGridLeafIterator<codim, All_Partition, const ThisType> (*this, true);
  }

  //! Iterator to first entity of codim 0 on leaf level (All_Partition)
  template <int dim, int dimworld>
  inline typename ALU2dGrid<dim, dimworld>::LeafIteratorType ALU2dGrid<dim, dimworld> :: leafbegin () const {
    return LeafIteratorImp(*this,false);
  }

  //! one past the end on this leaf level (codim 0 and All_Partition)
  template <int dim, int dimworld>
  inline typename ALU2dGrid<dim, dimworld>::LeafIteratorType ALU2dGrid<dim, dimworld>:: leafend () const {
    return LeafIteratorImp(*this,true);
  }

  //! for type identification
  template <int dim, int dimworld>
  inline GridIdentifier ALU2dGrid<dim, dimworld> :: type  () const {
    //type ALU2dGrid_Id has to be set in dune/common/grid.hh!
    return ALU3dGrid_Id;
  }

  //! Return maximum level defined in this grid. Levels are numbered
  //! 0 ... maxLevel with 0 the coarsest level.
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: maxLevel() const {
    return maxLevel_;
  }

  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::calcMaxlevel() {
    maxLevel_ = 0;
    ALU2DSPACE Listwalkptr <ALU2DSPACE Hmesh_basic::helement_t > walk(mesh_);
    for( walk->first() ; ! walk->done() ; walk->next()) {
      //Element & tr = walk->getitem();
      if(walk->getitem().level() > maxLevel_ )
        maxLevel_ = walk->getitem().level();
    }
  }


  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::calcExtras()
  {

    // Segmentation faults (!) beim SizeCache
    if(sizeCache_) delete sizeCache_;
    bool isSimplex = true;

    sizeCache_ = new SizeCacheType (*this,isSimplex,!isSimplex,true);

    for(unsigned int i=0; i<levelIndexVec_.size(); i++)
    {
      if(levelIndexVec_[i]) (*(levelIndexVec_[i])).calcNewIndex();
    }

    if(leafIndexSet_) leafIndexSet_->calcNewIndex();
    // update id set, i.e. insert new elements
    if(globalIdSet_) globalIdSet_->updateIdSet();
    //for(unsigned int i=0; i<MAXL; i++) vertexList_[i].unsetUp2Date();
    coarsenMarked_ = 0;
    refineMarked_  = 0;
  }

  //! Every time the grid is refined, data should be updated
  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::updateStatus() {
    calcMaxlevel();
    calcExtras();
  }

  //! number of grid entities in the entire grid for given codim
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld>:: hierSetSize (int cd) const {
    return mesh_.indexManagerSize(cd);
    //return size(cd);
  }

  //! number of grid entities per level and codim
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: size (int level, int cd) const
  {
    assert( sizeCache_ );
    return sizeCache_->size(level,cd);
  }

  //! number of entities per level, codim and geometry type in this process
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: size (int level, GeometryType type) const
  {
    if (type.isSimplex()) {
      return size(level, dim-type.dim());
    }
    return 0;
  }

  //! number of leaf entities per codim and geometry type in this process
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: size (GeometryType type) const
  {
    if (type.isSimplex()) {
      return size(dim-type.dim());
    }
    return 0;
  }

  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: size (int codim) const
  {
    assert( sizeCache_ );
    return sizeCache_->size(codim);
  }

  //! refine grid refCount times
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: globalRefine(int refCount) {

    for (int j = 0; j < refCount; ++j) {
      ALU2DSPACE Listwalkptr <ALU2DSPACE Hmesh_basic::helement_t > walk(mesh_);
      for( walk->first() ; ! walk->done() ; walk->next()) {
        ALU2DSPACE Element & tr = walk->getitem();
        tr.Refco_el::mark(ALU2DSPACE Refco::ref);
      }
      mesh_.refine();
    }
    postAdapt();
    //update data
    updateStatus();
    //hdl.refine() ist void!!!
    return true;
  }

  //! returns if a least one entity was marked for adaption
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: preAdapt () {
    return (refineMarked_  || coarsenMarked_);
  }

  //! clear all entity new markers
  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld> :: postAdapt () {
    ALU2DSPACE Listwalkptr <ALU2DSPACE Hmesh_basic::helement_t > walk(mesh_);
    for( walk->first() ; ! walk->done() ; walk->next()) {
      ALU2DSPACE Element & tr = walk->getitem();
      tr.ALU2DSPACE Refco_el::clear(ALU2DSPACE Refco::ref);
      tr.ALU2DSPACE Refco_el::clear(ALU2DSPACE Refco::crs);
    }
  }

  /**! refine all positive marked leaf entities,
     return true if a least one entity was refined
   */
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: adapt ( ) {
    if (preAdapt()) {
      mesh_.coarse();
      mesh_.refine();
      updateStatus();
      postAdapt();
      return true;
    }
    postAdapt();
    updateStatus();
    return false;
  }

  //! refine grid
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: refineGrid() {
    return adapt();
  }

  //! mark entities for refinement or coarsening, refCount < 0 will mark
  //! the entity for one coarsen step and refCount > 0 will mark for one
  //! refinement, one refinement will create 8 children per element
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: mark( int refCount , const typename Traits::template Codim<0>::EntityPointer & ep ) {
    return this->mark(refCount,*ep);
  }

  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: mark( int refCount , const typename Traits::template Codim<0>::Entity & en ) {
    //bool marked = (this->getRealImplementation(ep)).mark(ref);
    bool marked = true;
    if (refCount > 0)
      this->getRealImplementation(en).mark(ALU2DSPACE Refco::ref);
    else if (refCount < 0)
      this->getRealImplementation(en).mark(ALU2DSPACE Refco::crs);
    else marked = false;

    if(marked)
    {
      if(refCount > 0) refineMarked_ ++ ;
      if(refCount < 0) coarsenMarked_ ++ ;
    }
    return marked;
  }

  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::makeGeomTypes() {
    // stored is the dim, where is the codim
    for(int i=dim; i>= 0; i--)
      geomTypes_[dim-i][0] = GeometryType(GeometryType::simplex,i);
    return;
  }

  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::Traits :: LeafIndexSet &
  ALU2dGrid<dim, dimworld>::leafIndexSet() const
  {
    if(!leafIndexSet_) leafIndexSet_ = new LeafIndexSetImp ( *this );
    return *leafIndexSet_;
  }

  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::Traits :: LevelIndexSet &
  ALU2dGrid<dim, dimworld>::levelIndexSet( int level ) const
  {
    // check if level fits in vector
    assert( level >= 0 );
    assert( level < (int) levelIndexVec_.size() );

    if( levelIndexVec_[level] == 0 )
      levelIndexVec_[level] = new LevelIndexSetImp ( *this , level );
    return *(levelIndexVec_[level]);
  }


  template <int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld> & ALU2dGrid<dim, dimworld>::operator = (const ALU2dGrid<dim, dimworld> & g) {
    DUNE_THROW(GridError,"Do not use assignment operator of ALU2dGrid! \n");
    return (*this);
  }

  template <int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::ALU2dGrid(const ALU2dGrid<dim, dimworld> & g)
    : mesh_ (0)
      , maxLevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1,1)
      , hIndexSet_(*this)
      , globalIdSet_ (0)
      , localIdSet_ (*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
  {
    DUNE_THROW(GridError,"Do not use copy constructor of ALU2dGrid! \n");
  }

  template <int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::~ALU2dGrid()
  {
    for(unsigned int i=0; i<levelIndexVec_.size(); i++) delete levelIndexVec_[i];
    delete globalIdSet_; globalIdSet_ = 0;
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;
  }

}

#endif
