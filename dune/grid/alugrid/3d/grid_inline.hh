// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Dune includes
#include <dune/common/stdstreams.hh>

// Local includes
#include "entity.hh"
#include "iterator.hh"
#include "datahandle.hh"

namespace Dune {

  //--Grid
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3dGrid<dim, dimworld, elType>::
  ALU3dGrid(const std::string macroTriangFilename,
            const MPICommunicatorType mpiComm,
            const DuneBoundaryProjectionType* bndPrj,
            const DuneBoundaryProjectionVector* bndVec
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
      , referenceElement_( (elType == tetra) ?
                           GenericReferenceElements< alu3d_ctype, dim > :: simplex() :
                           GenericReferenceElements< alu3d_ctype, dim > :: cube() )
      , sizeCache_ (0)
      , lockPostAdapt_(false)
      , bndPrj_ ( bndPrj )
      , bndVec_ ( (bndVec) ? (new DuneBoundaryProjectionVector( *bndVec )) : 0 )
      , vertexProjection_( (bndPrj || bndVec) ? new ALUGridBoundaryProjectionType( *this ) : 0 )
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
      mygrid_ = createALUGrid( macroTriangFilename.c_str() );
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
  inline int ALU3dGrid<dim, dimworld, elType>::global_size(int codim) const
  {
    // return actual size of hierarchical index set
    // this is always up to date
    // maxIndex is the largest index used + 1
    return myGrid().indexManager(codim).getMaxIndex();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::hierSetSize(int codim) const
  {
    // return actual size of hierarchical index set
    return myGrid().indexManager(codim).getMaxIndex();
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline int ALU3dGrid<dim, dimworld, elType>::maxLevel() const
  {
    return maxlevel_;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline ALU3DSPACE GitterImplType* ALU3dGrid<dim, dimworld, elType>::
  createALUGrid(const char* macroName)
  {
    return new ALU3DSPACE GitterImplType (macroName
#if ALU3DGRID_PARALLEL
                                          , mpAccess_
#endif
#ifdef ALUGRID_VERTEX_PROJECTION
                                          , vertexProjection_ // only available in ALUGrid-1.15 or newer
#endif
                                          );
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
  inline typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::
  leafbegin(int level) const {
    return createLeafIteratorBegin<0, All_Partition> (level) ;
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
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
  inline typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
  ALU3dGrid<dim, dimworld, elType>::leafend(int level) const {
    return createLeafIteratorEnd <0, All_Partition> (level);
  }

  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline typename ALU3dGrid<dim, dimworld, elType>::LeafIteratorType
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

  // return Grid name
  template <int dim, int dimworld, ALU3dGridElementType elType>
  inline std::string ALU3dGrid<dim, dimworld, elType>::name ()
  {
    if(elType == hexa)
      return "ALUCubeGrid";
    assert( elType == tetra );
    return "ALUSimplexGrid";
  }

} // end namespace Dune
