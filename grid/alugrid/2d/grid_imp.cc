// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALU2DGRID_IMP_CC
#define DUNE_ALU2DGRID_IMP_CC

namespace Dune {

  //--Grid

  //! Constructor which reads an ALU2dGrid Macro Triang file
  //! or given GridFile
  template<int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::ALU2dGrid(std::string macroTriangFilename )
    :
#if ALU2DGRID_PARALLEL
      comm_( MPIHelper::getCommunicator() ),
#endif
      mygrid_ (new ALU2DSPACE Hmesh(checkMacroGridFile(macroTriangFilename)))
      , hIndexSet_(*this)
      , localIdSet_(*this)
      , levelIndexVec_(MAXL,0)
      , geomTypes_( dim+1 )
      , leafIndexSet_(0)
      , maxLevel_(0)
      , refineMarked_ (0)
      , coarsenMarked_ (0)
      , nrOfHangingNodes_(0)
      , sizeCache_(0)
      , lockPostAdapt_(false)
#if ALU2DGRID_PARALLEL
      , rankManager_( *this )
#endif
  {
#if ALU2DGRID_PARALLEL
    abort();
    rankManager_.initialize();
    //rankManager_.loadBalance();
#endif

    assert(mygrid_);
    makeGeomTypes();
    updateStatus();
  }

  template<int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::ALU2dGrid(std::string macroTriangFilename, int nrOfHangingNodes )
    :
#if ALU2DGRID_PARALLEL
      comm_( MPIHelper::getCommunicator() ),
#endif
      mygrid_ (new ALU2DSPACE Hmesh(checkMacroGridFile(macroTriangFilename),
                                    nrOfHangingNodes, ALU2DSPACE Refco::quart))
      , hIndexSet_(*this)
      , localIdSet_(*this)
      , levelIndexVec_( MAXL, 0 )
      , geomTypes_( dim+1 )
      , leafIndexSet_(0)
      , maxLevel_(0)
      , refineMarked_ (0)
      , coarsenMarked_ (0)
      , nrOfHangingNodes_(nrOfHangingNodes)
      , sizeCache_(0)
      , lockPostAdapt_(false)
#if ALU2DGRID_PARALLEL
      , rankManager_( *this )
#endif
  {
#if ALU2DGRID_PARALLEL
    rankManager_.initialize();
    //rankManager_.loadBalance();
#endif

    assert(mygrid_);
    makeGeomTypes();
    updateStatus();
  }

  //! Constructor which constructs an empty ALU2dGrid
  template<int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::ALU2dGrid(int nrOfHangingNodes)
    :
#if ALU2DGRID_PARALLEL
      comm_( MPIHelper::getCommunicator() ),
#endif
      mygrid_ (0)
      , hIndexSet_(*this)
      , localIdSet_(*this)
      , levelIndexVec_(MAXL,0)
      , geomTypes_( dim+1 )
      , leafIndexSet_(0)
      , maxLevel_(0)
      , refineMarked_ (0)
      , coarsenMarked_ (0)
      , nrOfHangingNodes_(nrOfHangingNodes)
      , sizeCache_(0)
      , lockPostAdapt_(false)
#if ALU2DGRID_PARALLEL
      , rankManager_( *this )
#endif
  {
    makeGeomTypes();
  }

  //! Iterator to first entity of given codim on level
  template <int dim, int dimworld>
  inline const char *
  ALU2dGrid<dim, dimworld> :: checkMacroGridFile(const std::string & filename)
  {
    std::ifstream file(filename.c_str());
    if(!file)
    {
      std::cerr << "Couldn't open file '" << filename <<"' !" << std::endl;
      DUNE_THROW(IOError,"Couldn't open file '" << filename <<"' !");
    }

    const std::string aluid("!Triangles");
    std::string idline;
    std::getline(file,idline);
    std::stringstream idstream(idline);
    std::string id;
    idstream >> id;

    if(id != aluid )
    {
      derr << "\nNon or wrong keyword '" << id << "' found! Expected keyword to be '" << aluid << "'! \n";
      DUNE_THROW(IOError,"Wrong file format!");
      return 0;
    }
    return filename.c_str();
  }



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

  //! for grid identification
  template <int dim, int dimworld>
  inline std::string ALU2dGrid<dim, dimworld> :: name () const
  {
    return ( nrOfHangingNodes_ > 0 ) ? "ALUSimplexGrid" : "ALUConformGrid";
  }

  //! Return maximum level defined in this grid. Levels are numbered
  //! 0 ... maxLevel with 0 the coarsest level.
  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: maxLevel() const {
    return maxLevel_;
  }

  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::calcMaxlevel()
  {
    maxLevel_ = 0;
    // walk the leaf level and take maximum as maxLevel
    ALU2DSPACE Listwalkptr <ALU2DSPACE Hmesh_basic::helement_t > walk( mesh() );
    for( walk->first() ; ! walk->done() ; walk->next())
    {
      if(walk->getitem().level() > maxLevel_ )
        maxLevel_ = walk->getitem().level();
    }
  }

  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::calcExtras()
  {
    // only in truly parallel runs update this
    //if( comm_.size() > 1 ) rankManager_.update();
#if ALU2DGRID_PARALLEL
    rankManager_.update();
#endif

    ////////////////////////////////////
    // reset size cache
    ////////////////////////////////////
    if(sizeCache_) delete sizeCache_;
    bool isSimplex = true;

    sizeCache_ = new SizeCacheType (*this,isSimplex,!isSimplex,true);
    /////////////////////////////////////

    // place this before update of level index, because
    // marker vector is used by level index set
    for(int i=0; i<MAXL; ++i) marker_[i].unsetUp2Date();
    leafMarker_.unsetUp2Date();

    ///////////////////////////////////////////////
    // update existing index sets
    ///////////////////////////////////////////////
    const int levelSize = levelIndexVec_.size();
    for(int i=0; i<levelSize; ++i)
    {
      if(levelIndexVec_[i]) (*(levelIndexVec_[i])).calcNewIndex();
    }

    if(leafIndexSet_) leafIndexSet_->calcNewIndex();
    ////////////////////////////////////////////////
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
    return mesh().indexManagerSize(cd);
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
  inline bool ALU2dGrid<dim, dimworld> :: globalRefine(int refCount)
  {
    if( refCount <= 0 ) return false;

    for (int j = 0; j < refCount; ++j)
    {
#if ALU2DGRID_PARALLEL
      {
        typedef typename Traits :: template Codim<0> :: template Partition<Interior_Partition> :: LeafIterator LeafIterator;
        LeafIterator endit = this->template leafend<0,Interior_Partition> ();
        for(LeafIterator it = this->template leafbegin<0,Interior_Partition> ();
            it != endit; ++it )
        {
          this->mark( 1, *it );
        }
      }

      rankManager_.notifyMarking();
#else
      ALU2DSPACE Listwalkptr <ALU2DSPACE Hmesh_basic::helement_t > walk(mesh());
      for( walk->first() ; ! walk->done() ; walk->next())
      {
        ALU2DSPACE Element & tr = walk->getitem();
        tr.Refco_el::mark(ALU2DSPACE Refco::ref);
      }
#endif
      mesh().refine();

      // in parallel update rank information
#if ALU2DGRID_PARALLEL
      rankManager_.update();
#endif
    }
    //update data
    updateStatus();

    // cleanup markers
    postAdapt();
    return true;
  }

  //! returns true if a least one entity was marked for coarseing
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: preAdapt () {
    return (coarsenMarked_ > 0);
  }


  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld> ::
  hierarchicClear (ALUElementType* el)
  {
    // clear actual tag
    el->Refco_el::clear();
    // clear refined tag
    el->Refco_el::clearWas();
    // go to children
    for(ALUElementType* child = el->down(); child; child = child->next())
    {
      // clear marker for child
      hierarchicClear(child);
    }
  }


  //! clear all entity new markers
  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld> :: postAdapt ()
  {
    // clear refinement markers throughout the grid
    typedef ALU2DSPACE Macro < ALU2DSPACE Element > macro_t;

    // get macro element iterator
    ALU2DSPACE Listwalkptr <macro_t> walk(mesh());
    for( walk->first() ; ! walk->done() ; walk->next())
    {
      // get element pointer
      ALUElementType* el = walk->getitem().operator ->();
      // hierarchically clear all markers
      hierarchicClear(el);
    }

    // reset marked element counters
    coarsenMarked_ = 0;
    refineMarked_  = 0;

    // unmark postAdapt flag
    lockPostAdapt_ = false;
  }


  /**! refine all positive marked leaf entities,
     return true if a least one entity was refined
   */
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: adapt ( )
  {
    if( lockPostAdapt_ == true )
    {
      DUNE_THROW(InvalidStateException,"Make sure that postAdapt is called after adapt was called and returned true!");
    }

#if ALU2DGRID_PARALLEL
    // make marking of ghost elements
    // needs one communication
    rankManager_.notifyMarking();
#else
    if( (refineMarked_ > 0) || (coarsenMarked_ > 0) )
#endif
    {
      // refine only will be done if
      // at least one element was marked for refinement
      bool adapted = (refineMarked_) ? true : false;
      mesh().refine();
      mesh().coarse();

      updateStatus();

      // notify that postAdapt must be called
      lockPostAdapt_ = true;

      return adapted;
    }
#if ALU2DGRID_PARALLEL == 0
    else
    {
      return false;
    }
#endif
  }


  // --adapt
  template< int dim, int dimworld >
  template< class GridImp, class DataHandle >
  inline bool ALU2dGrid< dim, dimworld >
  ::adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle )
  {
    typedef AdaptDataHandleInterface< GridImp, DataHandle > AdaptDataHandle;

    typedef typename EntityObject::ImplementationType EntityImp;
    EntityObject father( EntityImp(*this, this->maxLevel()) );
    EntityObject son( EntityImp(*this, this->maxLevel()) );

    int defaultChunk = newElementsChunk_;
    int actChunk     = refineEstimate_ * refineMarked_;

    // guess how many new elements we get
    int newElements = std::max( actChunk , defaultChunk );

    // reserve memory
    handle.preAdapt( newElements );

#if ALU2DGRID_PARALLEL
    // make marking of ghost elements
    // needs one communication
    rankManager_.notifyMarking();
#endif

    bool ref = false;
    {
      ALU2DSPACE AdaptRestrictProlong2dImpl< ThisType, AdaptDataHandle >
      rp(*this,
         father,this->getRealImplementation(father),
         son,   this->getRealImplementation(son),
         handle );

      ref = myGrid().duneAdapt(rp); // adapt grid
      //if(rp.maxLevel() >= 0) maxlevel_ = rp.maxLevel();
    }

    if(ref) {
      updateStatus();
    }

    // check whether we have balance
    handle.postAdapt();

    postAdapt();
    return ref;
  }


  // --adapt
  template <int dim, int dimworld>
  template <class DofManagerType, class RestrictProlongOperatorType>
  inline bool ALU2dGrid<dim, dimworld>::
  adapt(DofManagerType & dm, RestrictProlongOperatorType & rpo, bool verbose )
  {
    assert( ((verbose) ? (dverb << "ALU3dGrid :: adapt() new method called!\n", 1) : 1 ) );

    typedef RestrictProlongWrapper< ThisType, DofManagerType, RestrictProlongOperatorType > Wrapper;
    Wrapper rpOpWrapper( dm, rpo );
    const bool refined = adapt( rpOpWrapper );

    assert( ((verbose) ? (dverb << "ALU3dGrid :: adapt() new method finished!\n", 1) : 1 ) );
    return refined;
  }

  //! refine grid
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: refineGrid() {
    return adapt();
  }

  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> ::
  mark( int refCount , const typename Traits::template Codim<0>::Entity & en )
  {
    bool marked = this->getRealImplementation(en).mark(refCount);
    if(marked)
    {
      if(refCount > 0) ++refineMarked_;
      if(refCount < 0) ++coarsenMarked_;
    }
    return marked;
  }

  template <int dim, int dimworld>
  inline int ALU2dGrid<dim, dimworld> :: getMark(const typename Traits::template Codim<0>::Entity & e ) const
  {
    return this->getRealImplementation(e).getMark();
  }

  template <int dim, int dimworld>
  inline void ALU2dGrid<dim, dimworld>::makeGeomTypes()
  {
    // stored is the dim, where is the codim
    for(int i=dim; i>= 0; i--)
    {
      geomTypes_[ dim-i ].resize( 1 );
      geomTypes_[ dim-i ][ 0 ] = GeometryType( GeometryType :: simplex, i );
    }
  }

  //! get global id set of grid
  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::GlobalIdSet & ALU2dGrid<dim, dimworld>:: globalIdSet () const {
    return localIdSet();
  }

  //! get global id set of grid
  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::LocalIdSet & ALU2dGrid<dim, dimworld>:: localIdSet () const {
    return localIdSet_;
  }

  //! get hierarchic index set of the grid
  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::HierarchicIndexSet & ALU2dGrid<dim, dimworld>::hierarchicIndexSet () const {
    return hIndexSet_;
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

  // private methods that return underlying ALU2D Grid
  template <int dim, int dimworld>
  inline ALU2DSPACE Hmesh & ALU2dGrid<dim, dimworld>::myGrid()
  {
    return mesh();
  }
  template <int dim, int dimworld>
  inline ALU2DSPACE Hmesh & ALU2dGrid<dim, dimworld>::myGrid() const
  {
    return mesh();
  }

  //! return dummy communication
  template <int dim, int dimworld>
  inline const typename ALU2dGrid<dim, dimworld>::CollectiveCommunicationType & ALU2dGrid<dim, dimworld>::comm() const
  {
    return comm_;
  }


  template <int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::ALU2dGrid(const ALU2dGrid<dim, dimworld> & g)
    : mygrid_ (0)
      , maxLevel_(0)
      , coarsenMarked_(0) , refineMarked_(0)
      , geomTypes_(dim+1,1)
      , hIndexSet_(*this)
      , localIdSet_ (*this)
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
  {
    DUNE_THROW(GridError,"Do not use copy constructor of ALU2dGrid! \n");
  }

  template <int dim, int dimworld>
  inline ALU2dGrid<dim, dimworld>::~ALU2dGrid()
  {
    for(unsigned int i=0; i<levelIndexVec_.size(); i++) delete levelIndexVec_[i];
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;
    delete mygrid_;
  }

  // **************************************************************
  // ***************************************************************
  template <int dim, int dimworld>
  template <GrapeIOFileFormatType ftype>
  inline bool ALU2dGrid<dim, dimworld>::
  writeGrid(const std::string filename, alu2d_ctype time ) const
  {
    switch(ftype)
    {
    case xdr  : return writeGrid_Xdr(filename,time);
    case ascii : return writeGrid_Ascii(filename,time);
    default : derr << "Wrong file type in writeGrid method~ \n";
    }
    return false;
  }

  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld>::
  writeGrid_Ascii(const std::string filename, alu2d_ctype time ) const
  {
    abort();
    return true;
  }

  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld>::
  writeGrid_Xdr(const std::string filename, alu2d_ctype time ) const
  {
    ALU2DSPACE Hmesh & mygrd = myGrid();
    mygrd.storeGrid(filename.c_str(),time,0);

    // write time and maxlevel
    {
      std::string extraName(filename);
      extraName += ".extra";
      std::ofstream out (extraName.c_str());
      if(out)
      {
        out << std::scientific << time << " ";
        out << maxLevel_ << " ";
        out.close();
      }
      else
      {
        derr << "ALU2dGrid::writeGrid: couldn't open <" << extraName << ">! \n";
      }
    }
    return true;
  }

  template <int dim, int dimworld>
  template <GrapeIOFileFormatType ftype>
  inline bool ALU2dGrid<dim,dimworld>::
  readGrid( const std::string filename, alu2d_ctype & time )
  {
    {
      // if grid exists delete first
      if( mygrid_ ) delete mygrid_;
      mygrid_ = new ALU2DSPACE Hmesh (filename.c_str());
      assert(mygrid_ != 0);
    }
    {
      std::string extraName (filename);
      extraName += ".extra";
      std::ifstream in (extraName.c_str());
      if(in)
      {
        in  >> std::scientific >> time;
        in  >> maxLevel_;
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

    // cleanup markers
    postAdapt();
    return true;
  }


  // communicate level data
  template <int dim, int dimworld>
  template <class DataHandleImp,class DataType>
  inline void ALU2dGrid<dim, dimworld>::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir, int level ) const
  {
    // only communicate, if number of processes is larger than 1
#if ALU2DGRID_PARALLEL
    if( comm_.size() > 1 )
    {
      rankManager_.communicate(data,iftype,dir,level);
    }
#endif
  }

  // communicate level data
  template <int dim, int dimworld>
  template <class DataHandleImp,class DataType>
  inline void ALU2dGrid<dim, dimworld>::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir) const
  {
#if ALU2DGRID_PARALLEL
    // only communicate, if number of processes is larger than 1
    if( comm_.size() > 1 )
    {
      rankManager_.communicate(data,iftype,dir);
    }
#endif
  }

  // re-balance load of grid
  template <int dim, int dimworld>
  inline bool ALU2dGrid<dim, dimworld> :: loadBalance ()
  {
#if ALU2DGRID_PARALLEL
    if( comm_.size()  <= 1 ) return false;
    bool changed = rankManager_.loadBalance();
    if( changed ) updateStatus();
    return changed;
#else
    return false;
#endif
  }

  // re-balance load of grid
  template <int dim, int dimworld>
  template <class DataHandleImp>
  inline bool ALU2dGrid<dim, dimworld> :: loadBalance (DataHandleImp& data)
  {
#if ALU2DGRID_PARALLEL
    if( comm_.size()  <= 1 ) return false;
    bool changed = rankManager_.loadBalance();
    if( changed ) updateStatus();
    return changed;
#else
    return false;
#endif
  }
}

#endif
