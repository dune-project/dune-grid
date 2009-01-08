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
#include "treeiterator.cc"
#include "intersection.cc"

namespace Dune
{

  // AlbertaGrid::SetLocalCoords
  // ---------------------------

#ifndef CALC_COORD
  template< int dim, int dimworld >
  class AlbertaGrid< dim, dimworld >::SetLocalCoords
  {
    typedef Alberta::DofAccess< dim, dim > DofAccess;

    CoordVectorPointer coords_;
    DofAccess dofAccess_;

  public:
    explicit SetLocalCoords ( const CoordVectorPointer &coords )
      : coords_( coords ),
        dofAccess_( coords.dofSpace() )
    {}

    void operator() ( const Alberta::ElementInfo< dim > &elementInfo ) const
    {
      Alberta::GlobalVector *array = (Alberta::GlobalVector *)coords_;
      for( int i = 0; i < DofAccess::numSubEntities; ++i )
      {
        const Alberta::GlobalVector &x = elementInfo.coordinate( i );
        Alberta::GlobalVector &y = array[ dofAccess_( elementInfo.el(), i ) ];
        for( int i = 0; i < dimworld; ++i )
          y[ i ] = x[ i ];
      }
    }
  };
#endif



  // AlbertaGrid
  // -----------

  template< int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::AlbertaGrid ()
    : mesh_(),
      maxlevel_( 0 ),
      wasChanged_( false ),
      vertexMarkerLeaf_( false ), // creates LeafMarkerVector
      nv_( dim+1 ),
      dof_( 0 ),
      hIndexSet_( *this ),
      idSet_( *this ),
      levelIndexVec_( MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_ ( 0 ),
      coarsenMarked_( 0 ),
      refineMarked_( 0 ),
      lockPostAdapt_( false )
  {
    // If this check fails, reconfigure with --with-alberta-world-dim=dimworld
    dune_static_assert( (dimworld == DIM_OF_WORLD),
                        "Template Parameter dimworld does not match "
                        "ALBERTA's DIM_OF_WORLD setting." );

#if DUNE_ALBERTA_VERSION < 0x200
    // ALBERTA 1.2 supports only one grid dimension
    // If this check fails, reconfigure with --with-alberta-dim=dim
    dune_static_assert( (dim == DIM),
                        "Template Parameter dim does not match "
                        "ALBERTA's DIM setting." );
#endif
  }


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::initGrid ()
  {
    typedef Alberta::FillFlags< dim > FillFlags;

    dofNumbering_.create( mesh_ );

#ifndef CALC_COORD
    coords_.create( dofNumbering_.dofSpace( dimension ), "Coordinates" );
    SetLocalCoords setLocalCoords( coords_ );
    mesh_.hierarchicTraverse( setLocalCoords, FillFlags::coords );
    ((ALBERTA DOF_REAL_D_VEC *)coords_)->refine_interpol = &AlbertHelp::refineCoordsAndRefineCallBack< dimension  >;
    ((ALBERTA DOF_REAL_D_VEC *)coords_)->coarse_restrict = &AlbertHelp::coarseCallBack;
#endif

    elNewCheck_.create( dofNumbering_.dofSpace( 0 ), "el_new_check" );
    assert( !elNewCheck_ == false );
    elNewCheck_.template setupInterpolation< ElNewCheckInterpolation >();
    elNewCheck_.initialize( 0 );

    hIndexSet_.template initEntityNumbers< 0 >( dofNumbering_, elNumbers_[ 0 ] );
    hIndexSet_.template initEntityNumbers< 1 >( dofNumbering_, elNumbers_[ 1 ] );
    if( dim == 3 )
      hIndexSet_.template initEntityNumbers< 2 >( dofNumbering_, elNumbers_[ 2 ] );

    wasChanged_ = true;

    macroVertices_.resize( getMesh()->n_vertices );

    LeafDataType::initLeafDataValues( mesh_, 0 );

    calcExtras();
  }

  template < int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >
  ::AlbertaGrid ( const std::string &macroGridFileName,
                  const std::string &gridName )
    : mesh_( 0 ),
      maxlevel_( 0 ),
      wasChanged_( false ),
      vertexMarkerLeaf_( false ), // creates LeafMarkerVector
      nv_( dim+1 ),
      dof_( 0 ),
      hIndexSet_( *this ),
      idSet_( *this ),
      levelIndexVec_( MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_( 0 ),
      coarsenMarked_( 0 ),
      refineMarked_( 0 ),
      lockPostAdapt_( false )
  {
    // If this check fails, reconfigure with --with-alberta-world-dim=dimworld
    dune_static_assert( (dimworld == DIM_OF_WORLD),
                        "Template Parameter dimworld does not match "
                        "ALBERTA's DIM_OF_WORLD setting." );

#if DUNE_ALBERTA_VERSION < 0x200
    // ALBERTA 1.2 supports only one grid dimension
    // If this check fails, reconfigure with --with-alberta-dim=dim
    dune_static_assert( (dim == DIM),
                        "Template Parameter dim does not match "
                        "ALBERTA's DIM setting." );
#endif

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

    //ALBERTA AlbertHelp::initIndexManager_elmem_cc( hIndexSet_.indexStack_ );

    if( makeNew )
    {
      mesh_.create( macroGridFileName, gridName );
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


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::removeMesh ()
  {
    for( size_t i = 0; i < levelIndexVec_.size(); ++i )
    {
      if( levelIndexVec_[ i ] != 0 )
        delete levelIndexVec_[ i ];
      levelIndexVec_[ i ] = 0;
    }

    if( leafIndexSet_ != 0 )
      delete leafIndexSet_;
    leafIndexSet_ = 0;

    // release dof vectors
    for( int i = 0; i < AlbertHelp::numOfElNumVec; ++i )
      elNumbers_[ i ].release();
    elNewCheck_.release();
#ifndef CALC_COORD
    coords_.release();
#endif
    dofNumbering_.release();

    if( sizeCache_ != 0 )
      delete sizeCache_;
    sizeCache_ = 0;

    mesh_.release();
  }


  template< int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >::~AlbertaGrid ()
  {
    removeMesh();
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition<pitype>::LevelIterator
  AlbertaGrid< dim, dimworld >::lbegin ( int level ) const
  {
    assert( level >= 0 );

    if( level > maxlevel_ )
      return lend< codim, pitype >( level );

    if( codim > 0 )
    {
      if( !(vertexMarkerLevel_[ level ].up2Date() ) )
        vertexMarkerLevel_[ level ].markNewVertices( *this, level );
    }
    return AlbertaGridLevelIterator< codim, pitype, const This >( *this, &vertexMarkerLevel_[ level ], level );
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition< pitype >::LevelIterator
  AlbertaGrid < dim, dimworld >::lend ( int level ) const
  {
    assert( level >= 0 );
    return AlbertaGridLevelIterator< codim, pitype, const This >( *this, level );
  }


  template< int dim, int dimworld >
  template< int codim >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::LevelIterator
  AlbertaGrid < dim, dimworld >::lbegin ( int level ) const
  {
    return lbegin< codim, All_Partition >( level );
  }


  template< int dim, int dimworld >
  template< int codim >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::LevelIterator
  AlbertaGrid< dim, dimworld >::lend ( int level ) const
  {
    return lend< codim, All_Partition >( level );
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition< pitype >::LeafIterator
  AlbertaGrid< dim, dimworld >::leafbegin () const
  {
    if( codim > 1 )
    {
      if( !vertexMarkerLeaf_.up2Date() )
        vertexMarkerLeaf_.markNewLeafVertices( *this );
    }
    return AlbertaGridLeafIterator< codim, pitype, const This >( *this, &vertexMarkerLeaf_, maxlevel_ );
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition< pitype >::LeafIterator
  AlbertaGrid< dim, dimworld >::leafend () const
  {
    return AlbertaGridLeafIterator< codim, pitype, const This >( *this, maxlevel_ );
  }


  template< int dim, int dimworld >
  template< int codim >
  inline typename AlbertaGrid<dim,dimworld>::Traits
  ::template Codim< codim >::LeafIterator
  AlbertaGrid< dim, dimworld >::leafbegin () const
  {
    return leafbegin< codim, All_Partition >();
  }


  template< int dim, int dimworld >
  template< int codim >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::LeafIterator
  AlbertaGrid < dim, dimworld >::leafend () const
  {
    return leafend< codim, All_Partition >();
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
    typedef GetNewEntity< const This, EntityProvider, dim, codim > Helper;
    return Helper::getNewEntity( *this, entityProvider_ );
  }

  template< int dim, int dimworld >
  template< int codim >
  inline void AlbertaGrid< dim, dimworld >
  ::freeEntity ( typename SelectEntityImp< codim, dim, const This >::EntityObject * en ) const
  {
    typedef GetNewEntity< const This, EntityProvider, dim, codim > Helper;
    Helper::freeEntity( entityProvider_, en );
  }

  //**************************************
  //  refine and coarsen methods
  //**************************************

  template< int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::globalRefine ( int refCount )
  {
    typedef typename Traits::template Codim< 0 >::LeafIterator LeafIterator;

    // only MAXL level allowed
    assert( (refCount + maxlevel_) < MAXL );

    assert( refCount >= 0 );
    for( int i = 0; i < refCount; ++i )
    {
      // mark all interior elements
      const LeafIterator endit = leafend< 0 >();
      for( LeafIterator it = leafbegin< 0 >(); it != endit; ++it )
        mark( 1, *it );

      adapt();
      postAdapt();
    }

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
    Alberta::abs( elNewCheck_ );

    return wasChanged_;
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::mark( int refCount, const typename Traits::template Codim< 0 >::Entity &e )
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

    mark = refCount;
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


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >::adapt ()
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
    ALBERTA AlbertHelp::initIndexManager_elmem_cc( hIndexSet_.indexStack_ );

    // set all values of elNewCheck positive which means old
    Alberta::abs( elNewCheck_ );

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
    typedef typename SelectEntityImp< 0, dim, const This >::EntityImp EntityImp;

    EntityObject father( EntityImp( *this ) );
    EntityObject son( EntityImp( *this ) );

    typedef typename DofManagerType :: IndexSetRestrictProlongType IndexSetRPType;
    typedef CombinedAdaptProlongRestrict < IndexSetRPType,RestrictProlongOperatorType > COType;
    COType tmprpop ( dm.indexSetRPop() , data );

    int defaultElChunk = defaultElementChunk_;
    int newElementChunk_ = std::max( defaultElChunk , (refineMarked_ * 4) );

    // reserve memory
    dm.reserveMemory( newElementChunk_ );

    Alberta::AdaptRestrictProlongHandler< This, COType >
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
    maxlevel_ = Alberta::maxAbs( elNewCheck_ );
    assert( (maxlevel_ >= 0) && (maxlevel_ < MAXL) );

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
    typedef Alberta::FillFlags< dim > FillFlags;

    //removeMesh();

    if( filename.size() <= 0 )
      return false;

    mesh_.read( filename, time );
    if( !mesh_ )
      DUNE_THROW( AlbertaIOError, "Could not read grid file: " << filename << "." );

    dofNumbering_.create( mesh_ );

    for(int i=0; i<AlbertHelp::numOfElNumVec; i++)
    {
      std::ostringstream namestream;
      namestream << filename << "_num_c" << i;
      elNumbers_[ i ].read( namestream.str(), mesh_ );
    }

    elNewCheck_ = Alberta::DofVectorPointer< int >( AlbertHelp::getDofNewCheck( elNumbers_[ 0 ].dofSpace(), "el_new_check" ) );
    elNewCheck_.template setupInterpolation< ElNewCheckInterpolation >();

#ifndef CALC_COORD
    coords_.create( dofNumbering_.dofSpace( dimension ), "Coordinates" );
    SetLocalCoords setLocalCoords( coords_ );
    mesh_.hierarchicTraverse( setLocalCoords, FillFlags::coords );
    ((ALBERTA DOF_REAL_D_VEC *)coords_)->refine_interpol = &AlbertHelp::refineCoordsAndRefineCallBack< dimension  >;
    ((ALBERTA DOF_REAL_D_VEC *)coords_)->coarse_restrict = &AlbertHelp::coarseCallBack;
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
      hIndexSet_.indexStack_[i].setMaxIndex(maxIdx);
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


  template< int dim, int dimworld >
  struct AlbertaGrid< dim, dimworld >::ElNewCheckInterpolation
  {
    static const int dimension = dim;
    static const int codimension = 0;

  private:
    typedef Alberta::DofVectorPointer< int > DofVectorPointer;
    typedef Alberta::DofAccess< dimension, codimension > DofAccess;

    DofVectorPointer dofVector_;
    DofAccess dofAccess_;

    explicit ElNewCheckInterpolation ( const DofVectorPointer &dofVector )
      : dofVector_( dofVector ),
        dofAccess_( dofVector.dofSpace() )
    {}

  public:
    void operator() ( const Alberta::Element *father )
    {
      int *array = (int *)dofVector_;
      const int fatherLevel = std::abs( array[ dofAccess_( father, 0 ) ] );
      for( int i = 0; i < 2; ++i )
      {
        const Alberta::Element *child = father->child[ i ];
        array[ dofAccess_( child, 0 ) ] = -(fatherLevel+1);
      }
    }

    static void interpolateVector ( const DofVectorPointer &dofVector,
                                    const Alberta::Patch &patch )
    {
      ElNewCheckInterpolation interpolation( dofVector );
      patch.forEach( interpolation );
    }
  };

} // namespace Dune

#undef ALBERTA_CHAR
#endif
