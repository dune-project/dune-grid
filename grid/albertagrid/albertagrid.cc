// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_CC
#define DUNE_ALBERTAGRID_CC

//************************************************************************
//
//  implementation of AlbertaGrid
//
//  namespace Dune
//
//************************************************************************
#include "datahandle.hh"

#include "geometry.cc"
#include "entity.cc"
#include "entitypointer.cc"
#include "hierarchiciterator.cc"
#include "treeiterator.cc"
#include "intersection.cc"

namespace Dune
{

  template< int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::AlbertaGrid()
    : mesh_(),
      comm_(),
      maxlevel_( 0 ),
      wasChanged_( false ),
      vertexMarkerLeaf_( false ), // creates LeafMarkerVector
      nv_( dim+1 ),
      dof_( 0 ),
      hIndexSet_( *this ),
      globalIdSet_( *this ),
      levelIndexVec_( MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_ ( 0 ),
      coarsenMarked_( 0 ),
      refineMarked_( 0 ),
      lockPostAdapt_( false )
  {
    // create vector with geom types
    makeGeomTypes();
  }


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::makeGeomTypes ()
  {
    for( int codim = 0; codim <= dim; ++codim )
    {
      const GeometryType type( GeometryType::simplex, dim - codim );
      geomTypes_[ codim ].push_back( type );
    }
  }


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::initGrid ()
  {
    for( int i = 0; i < AlbertHelp::numOfElNumVec; ++i )
    {
      elNumbers_[ i ] = Alberta::DofVectorPointer< int >( AlbertHelp::getElNumbers( i ) );
      AlbertHelp::elNumbers[ i ] = NULL;
    }

    elNewCheck_ = Alberta::DofVectorPointer< int >( AlbertHelp::getElNewCheck() );
    AlbertHelp::elNewCheck = NULL;

#ifndef CALC_COORD
    coords_ = Alberta::DofVectorPointer< Alberta::GlobalVector >( AlbertHelp::getCoordVec< dimworld >() );
    AlbertHelp::coordVec = NULL;
#endif

    calcExtras();

    wasChanged_ = true;

    macroVertices_.resize( getMesh()->n_vertices );

    LeafDataType::initLeafDataValues( mesh_, 0 );

    calcExtras();
  }

  template < int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >
  ::AlbertaGrid ( const std::string &macroGridFileName, const std::string &gridName )
    : mesh_( 0 ),
      comm_(),
      maxlevel_( 0 ),
      wasChanged_( false ),
      vertexMarkerLeaf_( false ), // creates LeafMarkerVector
      nv_( dim+1 ),
      dof_( 0 ),
      hIndexSet_( *this ),
      globalIdSet_( *this ),
      levelIndexVec_( MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_( 0 ),
      coarsenMarked_( 0 ),
      refineMarked_( 0 ),
      lockPostAdapt_( false )
  {
    // create vector with geom types
    makeGeomTypes();

    if(dimworld != DIM_OF_WORLD)
    {
      DUNE_THROW(AlbertaError,"DUNE wasn't configured for dimworld = " <<
                 dimworld << ". Reconfigure with the --with-world-dim="<<dimworld<<" option!");
    }

    if(dim      != DIM)
    {
      DUNE_THROW(AlbertaError,"DUNE wasn't configured for dim = " <<
                 dim << ". Reconfigure with the --with-problem-dim="<<dim<<" option!");
    }

    bool makeNew = true;
    {
      std::fstream file( macroGridFileName.c_str(), std::ios::in );
      if( !file )
        DUNE_THROW( AlbertaIOError, "Could not open macro grid file '" << macroGridFileName << "'." );

      std::basic_string <char> str,str1;
      file >> str1; str = str1.assign(str1,0,3);
      // With that Albert MacroTriang starts DIM or DIM_OF_WORLD
      if (str != "DIM") makeNew = false;
      file.close();
    }

    ALBERTA AlbertHelp::initIndexManager_elmem_cc(indexStack_);

    if( makeNew )
    {
      mesh_.create( macroGridFileName, gridName );
      AlbertHelp::initDofAdmin< dim >( mesh_ );
      initGrid();
      std::cout << typeName() << " created from macro grid file '"
                << macroGridFileName << "'." << std::endl;
    }
    else
    {
      derr<<"Couldn't read grid file '"<< macroGridFileName<<"' because it's not in ALBERTA macro triangulation format! \n";
      DUNE_THROW(NotImplemented,"Constructor reading backup file not implemented!");
    }
  }

  template < int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::
  AlbertaGrid(const AlbertaGrid<dim,dimworld> & copy )
  {
    DUNE_THROW(AlbertaError,"do not use grid copy constructor! ");
  }

  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::removeMesh()
  {
    for(unsigned int i=0; i<levelIndexVec_.size(); i++)
    {
      if(levelIndexVec_[i])
      {
        delete levelIndexVec_[i];
        levelIndexVec_[i] = 0;
      }

    }

    if(leafIndexSet_)
    {
      delete leafIndexSet_;
      leafIndexSet_ = 0;
    }

    // release dof vectors
    for( int i = 0; i < AlbertHelp::numOfElNumVec; ++i )
      elNumbers_[ i ].release();
    elNewCheck_.release();
#ifndef CALC_COORD
    coords_.release();
#endif

    if(sizeCache_) delete sizeCache_;sizeCache_ = 0;

#if DIM == 3
    if(mesh_)
    {
      // because of bug in Alberta 1.2 , here until bug fixed
      ALBERTA RC_LIST_EL * rclist = ALBERTA get_rc_list(mesh_);
      rclist = 0;
    }
#endif
    mesh_.release();
  }

  // Desctructor
  template < int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::~AlbertaGrid()
  {
    removeMesh();
  }

  template < int dim, int dimworld >
  template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator
  //AlbertaGrid < dim, dimworld >::lbegin (int level, int proc) const
  AlbertaGrid < dim, dimworld >::lbegin (int level) const
  {
    assert( level >= 0 );
    // if we dont have this level return empty iterator
    if(level > maxlevel_) return this->template lend<codim,pitype> (level);

    if( codim > 0 ) //(dim == codim) || ((dim == 3) && (codim == 2)) )
    {
      if( ! (vertexMarkerLevel_[level].up2Date() ) )
        vertexMarkerLevel_[level].markNewVertices(*this,level);
    }
    return AlbertaGridLevelIterator<codim,pitype,const MyType> (*this,&vertexMarkerLevel_[level],level,-1);
  }

  template < int dim, int dimworld > template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator
  //AlbertaGrid < dim, dimworld >::lend (int level, int proc ) const
  AlbertaGrid < dim, dimworld >::lend (int level) const
  {
    return AlbertaGridLevelIterator<codim,pitype,const MyType> ((*this),level,-1);
  }

  template < int dim, int dimworld > template<int codim>
  inline typename AlbertaGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
  AlbertaGrid < dim, dimworld >::lbegin (int level) const
  {
    return this->template lbegin<codim,All_Partition> (level);
  }

  template < int dim, int dimworld > template<int codim>
  inline typename AlbertaGrid<dim, dimworld>::Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
  AlbertaGrid < dim, dimworld >::lend (int level) const
  {
    return this->template lend<codim,All_Partition> (level);
  }

  template < int dim, int dimworld >
  template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin (int level, int proc ) const
  {
    if((dim == codim) || ((dim == 3) && (codim == 2)) )
    {
      if( ! (vertexMarkerLeaf_.up2Date()) ) vertexMarkerLeaf_.markNewLeafVertices(*this);
    }
    return AlbertaGridLeafIterator<codim, pitype, const MyType> (*this,&vertexMarkerLeaf_,level,proc);
  }

  template < int dim, int dimworld >
  template<int codim>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin (int level, int proc ) const {
    return leafbegin<codim, All_Partition>(level, proc);
  }


  template < int dim, int dimworld >
  template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin () const {
    return leafbegin<codim, pitype>(maxlevel_, -1);
  }

  template < int dim, int dimworld >
  template<int codim>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin () const {
    return leafbegin<codim, All_Partition>(maxlevel_, -1);
  }


  template < int dim, int dimworld >
  template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend (int level, int proc ) const
  {
    return AlbertaGridLeafIterator<codim, pitype, const MyType> (*this,level,proc);
  }

  template < int dim, int dimworld >
  template<int codim>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend (int level, int proc ) const {
    return leafend<codim, All_Partition>(level, proc);
  }

  template < int dim, int dimworld >
  template<int codim, PartitionIteratorType pitype>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend () const {
    return leafend<codim, pitype>(maxlevel_, -1);
  }

  template < int dim, int dimworld >
  template<int codim>
  inline typename AlbertaGrid<dim,dimworld>::Traits::template Codim<codim>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend () const {
    return leafend<codim, All_Partition>(maxlevel_, -1);
  }

  template < int dim, int dimworld >
  inline typename AlbertaGrid<dim,dimworld>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin (int level, int proc ) const
  {
    return leafbegin<0, All_Partition> (level,proc);
  }

  template < int dim, int dimworld >
  inline typename AlbertaGrid<dim,dimworld>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafbegin () const {
    return leafbegin<0, All_Partition>(maxlevel_, -1);
  }

  template < int dim, int dimworld >
  inline typename AlbertaGrid<dim,dimworld>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend (int level, int proc ) const
  {
    return leafend<0,All_Partition> (level,proc);
  }

  template < int dim, int dimworld >
  inline typename AlbertaGrid<dim,dimworld>::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend () const {
    return leafend<0, All_Partition>(maxlevel_, -1);
  }

  //****************************************
  // getNewEntity methods
  //***************************************

  // default implementation used new and delete
  template< class GridImp, class EntityProvider, int dim , int codim >
  struct GetNewEntity
  {
    typedef typename SelectEntityImp<codim,dim,GridImp>::EntityObject EntityObject;
    typedef typename SelectEntityImp<codim,dim,GridImp>::EntityImp EntityImp;

    static EntityObject *
    getNewEntity ( GridImp &grid, EntityProvider &enp )
    {
      return new EntityObject( EntityImp( grid ) );
    }

    static void freeEntity ( EntityProvider &enp , EntityObject *en )
    {
      if( en )
        delete en;
    }
  };

  // specialisation for codim 0 uses stack
  template <class GridImp, class EntityProvider, int dim>
  struct GetNewEntity<GridImp,EntityProvider,dim,0>
  {
    typedef typename SelectEntityImp<0,dim,GridImp>::EntityObject EntityObject;
    typedef typename SelectEntityImp<0,dim,GridImp>::EntityImp EntityImp;

    static EntityObject *
    getNewEntity ( GridImp &grid, EntityProvider &enp )
    {
      // return object from stack
      return enp.getNewObjectEntity( grid, (EntityImp *)0 );
    }

    static void freeEntity ( EntityProvider &enp , EntityObject *en )
    {
      enp.freeObjectEntity( en );
    }
  };

  template< int dim, int dimworld >
  template< int codim >
  inline typename
  SelectEntityImp< codim, dim, const AlbertaGrid< dim, dimworld > >::EntityObject *
  AlbertaGrid< dim, dimworld >::getNewEntity () const
  {
    typedef GetNewEntity< const MyType, EntityProvider, dim, codim > Helper;
    return Helper::getNewEntity( *this, entityProvider_ );
  }

  template< int dim, int dimworld >
  template< int codim >
  inline void AlbertaGrid< dim, dimworld >
  ::freeEntity ( typename SelectEntityImp< codim, dim, const MyType >::EntityObject * en ) const
  {
    typedef GetNewEntity< const MyType, EntityProvider, dim, codim > Helper;
    Helper::freeEntity( entityProvider_, en );
  }

  //**************************************
  //  refine and coarsen methods
  //**************************************
  //  --Adaptation
  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::
  globalRefine(int refCount)
  {
    // only MAXL level allowed
    assert( (refCount + maxlevel_) < MAXL );

    typedef LeafIterator LeafIt;
    LeafIt endit = this->leafend(this->maxLevel());

    assert(refCount >= 0);
    for(int i=0; i<refCount; i++)
    {
      // mark all interior elements
      for(LeafIt it = this->leafbegin(this->maxLevel()); it != endit; ++it)
      {
        this->mark(1, *it );
      }

      // mark all ghosts
      for(LeafIt it = leafbegin(this->maxLevel(),Ghost_Partition); it != endit; ++it)
      {
        this->mark(1, *it );
      }

      this->adapt();
      this->postAdapt();
    }

    //std::cout << "Grid was global refined !\n";
    return wasChanged_;
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >::preAdapt ()
  {
    return (coarsenMarked_ > 0);
  }


  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::postAdapt()
  {
    typedef Alberta::Mesh Mesh;
    assert( (leafIndexSet_) ? (((Mesh *)mesh_)->n_elements == leafIndexSet_->size(0) ?   1 : 0) : 1);
    assert( (leafIndexSet_) ? (((Mesh *)mesh_)->n_vertices == leafIndexSet_->size(dim) ? 1 : 0) : 1);
#if DIM == 3
    //assert( (leafIndexSet_ && dim == 3) ? (mesh->n_edges == leafIndexSet_->size(dim-1) ?  1 :0) :1);
    assert( (leafIndexSet_ && dim == 3) ? (((Mesh *)mesh_)->n_faces == leafIndexSet_->size(1) ? 1 : 0) : 1);
#endif
    // if lockPostAdapt == false, the user forgot to call adapt before postAdapt
    if( lockPostAdapt_ == false )
    {
      DUNE_THROW(InvalidStateException,"AlbertaGrid::postAdapt called without previous adapt call!");
    }

    // unlock post adapt
    lockPostAdapt_ = false;

    coarsenMarked_ = 0;
    refineMarked_  = 0;

    // clear refined marker
    ALBERTA AlbertHelp::set2positive( elNewCheck_ );

    return wasChanged_;
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::mark( int refCount, const typename Traits::template Codim< 0 >::Entity &e ) const
  {
    // if not leaf entity, leave method
    if( !e.isLeaf() )
      return false;

    // take back previous marking
    int mark = getMark( e );
    if( mark < 0 )
      --coarsenMarked_;
    if( mark > 0 )
      refineMarked_ -= (2 << mark);

    // set new marking (max( previous, refCount ))
    mark = std::max( mark, refCount );
    if( mark < 0 )
      ++coarsenMarked_;
    if( mark > 0 )
      refineMarked_ += (2 << mark);
    getRealImplementation( e ).elementInfo().setMark( mark );

    return true;
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >
  ::getMark( const typename Traits::template Codim< 0 >::Entity &e ) const
  {
    return getRealImplementation( e ).elementInfo().getMark();
  }

  // --adapt
  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::adapt()
  {
    wasChanged_ = false;

    // if lockPostAdapt == true, the user forgot to call postAdapt
    // in previous cycles
    if( lockPostAdapt_ == true )
    {
      DUNE_THROW(InvalidStateException,"AlbertaGrid::adapt called without previous postAdapt call!");
    }
    // lock for post Adapt
    lockPostAdapt_ = true;

    // set global pointer to index manager in elmem.cc
    ALBERTA AlbertHelp::initIndexManager_elmem_cc( indexStack_ );

    // set all values of elNewCheck positive which means old
    ALBERTA AlbertHelp::set2positive ( elNewCheck_ );

    const bool refined = mesh_.refine();
    const bool coarsened = (preAdapt() ? mesh_.coarsen() : false);
    wasChanged_ = (refined || coarsened);

    if( wasChanged_ )
      calcExtras();

    // remove global pointer in elmem.cc
    ALBERTA AlbertHelp::removeIndexManager_elmem_cc(AlbertHelp::numOfElNumVec);

    // return true if elements were created
    return refined;
  }

  template < int dim, int dimworld >
  template <class DofManagerType, class RestrictProlongOperatorType>
  inline bool AlbertaGrid < dim, dimworld >::
  adapt(DofManagerType & dm, RestrictProlongOperatorType & data, bool verbose)
  {
#ifndef CALC_COORD
    typedef typename SelectEntityImp<0,dim,const MyType>::EntityImp EntityImp;

    EntityObject father(EntityImp(*this,maxLevel(),true));
    EntityObject son(EntityImp(*this,maxLevel(),true));

    typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
    typedef CombinedAdaptProlongRestrict < IndexSetRPType,RestrictProlongOperatorType > COType;
    COType tmprpop ( dm.indexSetRPop() , data );

    int defaultElChunk = defaultElementChunk_;
    int newElementChunk_ = std::max( defaultElChunk , (refineMarked_ * 4) );

    // reserve memory
    dm.reserveMemory( newElementChunk_ );

    Alberta::AdaptRestrictProlongHandler< MyType , COType >
    handler ( *this, father, this->getRealImplementation( father ),
              son, this->getRealImplementation( son ), tmprpop );

    ALBERTA AlbertHelp::MeshCallBack &callBack = ALBERTA
                                                 AlbertHelp::MeshCallBack::instance();

    callBack.setPointers( mesh_, handler );

    bool refined = this->adapt();

    callBack.reset();

    dm.dofCompress();
    postAdapt();
    return refined;
#else
    derr << "AlbertaGrid::adapt(dm,rp) : CALC_COORD defined, therefore adapt with call back not defined! \n";
    return false ;
#endif
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::checkElNew ( const Alberta::Element *element ) const
  {
    // if element is new then entry in dofVec is 1
    const int *elNewVec = (int *)elNewCheck_;
    return (elNewVec[ element->dof[ dof_ ][ nv_ ] ] < 0);
  }


  template< int dim, int dimworld >
  inline const Alberta::GlobalVector &
  AlbertaGrid< dim, dimworld >
  ::getCoord ( const Alberta::ElementInfo< dimension > &elementInfo,
               int vertex ) const
  {
    assert( (vertex >= 0) && (vertex <= dim) );
#ifdef CALC_COORD
    return elementInfo.coordinate( vx );
#else
    assert( !(!coords_) );
    Alberta::GlobalVector *coordsVec = (Alberta::GlobalVector *)coords_;
    return coordsVec[ elementInfo.el()->dof[ vertex ][ 0 ] ];
#endif
  }


  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::maxLevel() const
  {
    return maxlevel_;
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::global_size (int codim) const
  {
    if(codim == dim) return getMesh()->n_vertices;
    // for higher codims we have the index stack
    return indexStack_[codim].size();
  }

  // --size
  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::size (int level, int codim) const
  {
    if( (level > maxlevel_) || (level < 0) ) return 0;
    assert( sizeCache_ );
    return sizeCache_->size(level,codim);
  }


  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::size (int level, GeometryType type) const
  {
    return ((type.isSimplex()) ? this->size(level,dim-type.dim()) : 0);
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::size (GeometryType type) const
  {
    return ((type.isSimplex()) ? this->size(dim-type.dim()) : 0);
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::size (int codim) const
  {
    assert( sizeCache_ );
    return sizeCache_->size(codim);
  }

  template < int dim, int dimworld >
  inline const typename AlbertaGrid < dim, dimworld > :: Traits :: LevelIndexSet &
  AlbertaGrid < dim, dimworld > :: levelIndexSet (int level) const
  {
    // assert that given level is in range
    assert( level >= 0 );
    assert( level < (int) levelIndexVec_.size() );

    if(!levelIndexVec_[level]) levelIndexVec_[level] = new LevelIndexSetImp (*this,level);
    return *(levelIndexVec_[level]);
  }

  template < int dim, int dimworld >
  inline const typename AlbertaGrid < dim, dimworld > :: Traits :: LeafIndexSet &
  AlbertaGrid < dim, dimworld > :: leafIndexSet () const
  {
    if(!leafIndexSet_) leafIndexSet_ = new LeafIndexSetImp (*this);
    return *leafIndexSet_;
  }


  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::arrangeDofVec()
  {
    hIndexSet_.updatePointers( elNumbers_ );

    elAdmin_ = elNumbers_[ 0 ].dofSpace()->admin;

    // see Albert Doc. , should stay the same
    const_cast<int &> (nv_)  = elAdmin_->n0_dof[CENTER];
    const_cast<int &> (dof_) = elAdmin_->mesh->node[CENTER];
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >
  ::getLevelOfElement ( const Alberta::Element *element ) const
  {
    assert( element != NULL );
    const int *elNewVec = (int *)elNewCheck_;
    // return the elements level which is the absolute value of the entry
    return std::abs( elNewVec[ element->dof[ dof_ ][ nv_ ] ] );
  }


  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::getElementNumber ( const ALBERTA EL * el ) const
  {
    return hIndexSet_.getIndex(el,0,Int2Type<dim>());
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::getFaceNumber ( const ALBERTA EL * el , int face ) const
  {
    // codim of faces is 2 therefore dim-1
    assert( face >= 0 );
    assert( face < dim+1 );
    return hIndexSet_.getIndex(el,face,Int2Type<dim-1>());
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::getEdgeNumber ( const ALBERTA EL * el , int edge ) const
  {
    assert(dim == 3);
    // codim of edges is 2 therefore dim-2
    return hIndexSet_.getIndex(el,edge,Int2Type<dim-2>());
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::getVertexNumber ( const ALBERTA EL * el , int vx ) const
  {
    return hIndexSet_.getIndex(el,vx,Int2Type<0>());
  }

  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::calcExtras ()
  {
    // store pointer to numbering vectors and so on
    arrangeDofVec ();

    // determine new maxlevel
    maxlevel_ = ALBERTA AlbertHelp::calcMaxAbsoluteValueOfVector( elNewCheck_ );
    assert( maxlevel_ >= 0);
    assert( maxlevel_ < MAXL);

#ifndef NDEBUG
    int mlvl = ALBERTA AlbertHelp::calcMaxLevel( mesh_, elNewCheck_ );
    assert( mlvl == maxlevel_ );
#endif

    // unset up2Dat status, if lbegin is called then this status is updated
    for(int l=0; l<MAXL; ++l) vertexMarkerLevel_[l].unsetUp2Date();

    // unset up2Dat status, if leafbegin is called then this status is updated
    vertexMarkerLeaf_.unsetUp2Date();

    if(sizeCache_) delete sizeCache_;
    // first bool says we have simplex, second not cube, third, worryabout
    sizeCache_ = new SizeCacheType (*this,true,false,true);

    // if levelIndexSet exists, then update now
    for(unsigned int i=0; i<levelIndexVec_.size(); ++i)
      if(levelIndexVec_[i]) (*levelIndexVec_[i]).calcNewIndex();

    // create new Leaf Index
    if( leafIndexSet_ ) leafIndexSet_->calcNewIndex();

    // we have a new grid
    wasChanged_ = true;
  }


  template< int dim, int dimworld >
  template< GrapeIOFileFormatType format >
  inline bool AlbertaGrid< dim, dimworld >
  ::writeGrid ( const std::string &filename, ctype time ) const
  {
    switch( format )
    {
    case xdr :
      return writeGridXdr( filename, time );

    case ascii :
      DUNE_THROW( NotImplemented, "AlbertaGrid does not support writeGrid< ascii >." );

    // write leaf grid as macro triangulation
    //int ret = ALBERTA write_macro( mesh_ , filename.c_str() );
    //return (ret == 1) ? true : false;

    default :
      DUNE_THROW( NotImplemented, "AlbertaGrid: Unknown output format: " << format << "." );
    }
    return false;
  }


  template< int dim, int dimworld >
  template< GrapeIOFileFormatType format >
  inline bool AlbertaGrid< dim, dimworld >
  ::readGrid ( const std::string &filename, ctype &time )
  {
    switch( format )
    {
    case xdr :
      return readGridXdr( filename, time );

    case ascii :
      DUNE_THROW( NotImplemented, "AlbertaGrid does not support readGrid< ascii >." );

    //return readGridAscii (filename , time );

    default :
      DUNE_THROW( NotImplemented, "AlbertaGrid: Unknown output format: " << format << "." );
    }
    return false;
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::writeGridXdr ( const std::string &filename, ctype time ) const
  {
    if( filename.size() <= 0 )
      DUNE_THROW( AlbertaIOError, "No filename given to writeGridXdr." );

    bool success = mesh_.write( filename, time );

    // strore element numbering to file
    for( int i = 0; i < AlbertHelp::numOfElNumVec; ++i )
    {
      std::ostringstream namestream;
      namestream << filename << "_num_c" << i;
      success &= elNumbers_[ i ].write( namestream.str() );
    }

    return success;
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::readGridXdr ( const std::string &filename, ctype &time )
  {
    // remove all old stuff
    // to be reivised
    //removeMesh();

    if( filename.size() <= 0 )
      return false;

    mesh_.read( filename, time );
    if( !mesh_ )
      DUNE_THROW( AlbertaIOError, "Could not read grid file: " << filename << "." );

    for(int i=0; i<AlbertHelp::numOfElNumVec; i++)
    {
      std::ostringstream namestream;
      namestream << filename << "_num_c" << i;
      elNumbers_[ i ].read( namestream.str(), mesh_ );
    }

    elNewCheck_ = Alberta::DofVectorPointer< int >( AlbertHelp::getDofNewCheck( elNumbers_[ 0 ].dofSpace(), "el_new_check" ) );

#ifndef CALC_COORD
    assert( !coords_ );
    coords_ = Alberta::DofVectorPointer< int >( AlbertHelp::makeTheRest< dim, dimworld >( mesh_ ) );
#endif

    // restore level information for each element by traversing the mesh
    ALBERTA AlbertHelp::restoreElNewCheck( mesh_ , elNewCheck_ );

    // make vectors know in grid and hSet
    arrangeDofVec();

    // calc maxlevel and indexOnLevel and so on
    calcExtras();

    // set el_index of index manager to max element index
    for(int i=0; i<ALBERTA AlbertHelp::numOfElNumVec; i++)
    {
      int maxIdx = ALBERTA AlbertHelp::calcMaxIndex( elNumbers_[ i ] );
      indexStack_[i].setMaxIndex(maxIdx);
    }

    return true;
  }


#if 0
  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::readGridAscii
    (const std::basic_string<char> filename, albertCtype & time )
  {
    removeMesh(); // delete all objects

    mesh_.create( "AlbertaGrid", filename.c_str() );
    AlbertHelp::initDofAdmin< dim >( mesh_ );

    time = 0.0;

    // unset up2Dat status, if lbegin is called then this status is updated
    for(int l=0; l<MAXL; l++) vertexMarkerLevel_[l].unsetUp2Date();

    // unset up2Dat status, if leafbegin is called then this status is updated
    vertexMarkerLeaf_.unsetUp2Date();

    ALBERTA AlbertHelp::initIndexManager_elmem_cc(indexStack_);

    initGrid();
    return true;
  }
#endif

} // namespace Dune

#undef ALBERTA_CHAR
#endif
