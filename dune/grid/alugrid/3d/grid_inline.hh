// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Dune includes
#include <dune/common/stdstreams.hh>

// Local includes
#include "alu3dinclude.hh"
#include "entity.hh"
#include "iterator.hh"
#include "datahandle.hh"
#include "grid.hh"

namespace Dune
{

  // Implementation of ALU3dGrid
  // ---------------------------

  template< ALU3dGridElementType elType, class Comm >
  inline ALU3dGrid< elType, Comm >
  ::ALU3dGrid ( const std::string &macroTriangFilename,
                const MPICommunicatorType mpiComm,
                const DuneBoundaryProjectionType *bndPrj,
                const DuneBoundaryProjectionVector *bndVec,
                const ALUGridRefinementType refinementType )
    : mygrid_( 0 )
      , maxlevel_( 0 )
      , coarsenMarked_( 0 )
      , refineMarked_( 0 )
      , geomTypes_() //dim+1, std::vector<GeometryType>(1) )
      , hIndexSet_ (*this)
      , globalIdSet_( 0 )
      , localIdSet_( *this )
      , levelIndexVec_(MAXL,0) , leafIndexSet_(0)
      , referenceElement_( elType == tetra
                           ? GenericReferenceElements< alu3d_ctype, dimension > :: simplex()
                           : GenericReferenceElements< alu3d_ctype, dimension > :: cube() )
      , sizeCache_ ( 0 )
#ifdef USE_SMP_PARALLEL
      , factoryVec_( GridObjectFactoryType :: maxThreads(), GridObjectFactoryType( *this ) )
#else
      , factory_( *this )
#endif
      , lockPostAdapt_( false )
      , bndPrj_ ( bndPrj )
      , bndVec_ ( (bndVec) ? (new DuneBoundaryProjectionVector( *bndVec )) : 0 )
      , vertexProjection_( (bndPrj || bndVec) ? new ALUGridBoundaryProjectionType( *this ) : 0 )
      , communications_( new Communications( mpiComm ) )
      , refinementType_( refinementType )
  {
    assert( elType == tetra || elType == hexa );

    geomTypes_.resize( dimension+1 );
    GeometryType tmpType;
    for( int codim = 0; codim <= dimension; ++codim )
    {
      if (elType == tetra)
        tmpType.makeSimplex( dimension - codim );
      else
        tmpType.makeCube( dimension - codim );

      geomTypes_[ codim ].push_back( tmpType );
    }

    // check macro grid file for keyword
    checkMacroGridFile( macroTriangFilename );

    mygrid_ = createALUGrid( macroTriangFilename );
    assert( mygrid_ );

    dverb << "************************************************" << std::endl;
    dverb << "Created grid on p=" << comm().rank() << std::endl;
    dverb << "************************************************" << std::endl;
    checkMacroGrid ();

    postAdapt();
    calcExtras();
  } // end constructor


  template< ALU3dGridElementType elType, class Comm >
  inline int ALU3dGrid< elType, Comm >::global_size ( int codim ) const
  {
    // return actual size of hierarchical index set
    // this is always up to date
    // maxIndex is the largest index used + 1
    return myGrid().indexManager(codim).getMaxIndex();
  }


  template< ALU3dGridElementType elType, class Comm >
  inline int ALU3dGrid< elType, Comm >::hierSetSize ( int codim ) const
  {
    // return actual size of hierarchical index set
    return myGrid().indexManager(codim).getMaxIndex();
  }


  template< ALU3dGridElementType elType, class Comm >
  inline int ALU3dGrid< elType, Comm >::maxLevel () const
  {
    return maxlevel_;
  }


  template< ALU3dGridElementType elType, class Comm >
  inline typename ALU3dGrid< elType, Comm >::GitterImplType &
  ALU3dGrid< elType, Comm >::myGrid () const
  {
    assert( mygrid_ );
    return *mygrid_;
  }


  // lbegin methods
  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LevelIterator
  ALU3dGrid< elType, Comm >::lbegin ( int level ) const
  {
    assert( level >= 0 );
    // if we dont have this level return empty iterator
    if( level > maxlevel_ )
      return this->template lend<cd,pitype> (level);

    return ALU3dGridLevelIterator< cd, pitype, const ThisType >( factory(), level, true );
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LevelIterator
  ALU3dGrid< elType, Comm >::lend ( int level ) const
  {
    assert( level >= 0 );
    return ALU3dGridLevelIterator< cd, pitype, const ThisType >( factory(), level );
  }


  // lbegin methods
  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< All_Partition >::LevelIterator
  ALU3dGrid< elType, Comm >::lbegin ( int level ) const
  {
    return this->template lbegin<cd,All_Partition>( level );
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< All_Partition >::LevelIterator
  ALU3dGrid< elType, Comm >::lend ( int level ) const
  {
    assert( level >= 0 );
    return this->template lend<cd,All_Partition>( level );
  }


  //***********************************************************
  //
  // leaf methods , first all begin methods
  //
  //***********************************************************
  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::createLeafIteratorBegin ( int level ) const
  {
    assert( level >= 0 );
    return ALU3dGridLeafIterator< cd, pitype, const ThisType >( factory(), level, true );
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::leafbegin ( int level ) const
  {
    return createLeafIteratorBegin<cd, pitype> (level) ;
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::LeafIterator
  ALU3dGrid< elType, Comm >::leafbegin ( int level ) const
  {
    return createLeafIteratorBegin<cd, All_Partition> (level) ;
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::leafbegin () const
  {
    return createLeafIteratorBegin< cd, pitype > (maxlevel_) ;
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::LeafIterator
  ALU3dGrid< elType, Comm >::leafbegin () const
  {
    return createLeafIteratorBegin< cd, All_Partition> (maxlevel_) ;
  }


  template< ALU3dGridElementType elType, class Comm >
  inline typename ALU3dGrid< elType, Comm >::LeafIteratorType
  ALU3dGrid< elType, Comm >::leafbegin ( int level ) const
  {
    return createLeafIteratorBegin<0, All_Partition> (level) ;
  }


  template< ALU3dGridElementType elType, class Comm >
  inline typename ALU3dGrid< elType, Comm >::LeafIteratorType
  ALU3dGrid< elType, Comm >::leafbegin () const
  {
    return createLeafIteratorBegin<0, All_Partition> (maxlevel_) ;
  }


  //****************************************************************
  //
  // all leaf end methods
  //
  //****************************************************************
  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::createLeafIteratorEnd ( int level ) const
  {
    assert( level >= 0 );
    return ALU3dGridLeafIterator<cd, pitype, const MyType> ( factory() , level);
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::leafend ( int level ) const
  {
    return createLeafIteratorEnd < cd, pitype> (level);
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::LeafIterator
  ALU3dGrid< elType, Comm >::leafend ( int level ) const
  {
    return createLeafIteratorEnd < cd, All_Partition> (level);
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd, PartitionIteratorType pitype >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::template Partition< pitype >::LeafIterator
  ALU3dGrid< elType, Comm >::leafend () const
  {
    return createLeafIteratorEnd < cd, pitype> (maxlevel_);
  }


  template< ALU3dGridElementType elType, class Comm >
  template< int cd >
  inline typename ALU3dGrid< elType, Comm >::Traits::template Codim< cd >::LeafIterator
  ALU3dGrid< elType, Comm >::leafend () const
  {
    return createLeafIteratorEnd < cd, All_Partition> (maxlevel_);
  }


  template< ALU3dGridElementType elType, class Comm >
  inline typename ALU3dGrid< elType, Comm >::LeafIteratorType
  ALU3dGrid< elType, Comm >::leafend ( int level ) const
  {
    return createLeafIteratorEnd <0, All_Partition> (level);
  }


  template< ALU3dGridElementType elType, class Comm >
  inline typename ALU3dGrid< elType, Comm >::LeafIteratorType
  ALU3dGrid< elType, Comm >::leafend () const
  {
    return createLeafIteratorEnd <0,All_Partition> (maxlevel_);
  }


  //*****************************************************************

  // mark given entity
  template< ALU3dGridElementType elType, class Comm >
  inline bool ALU3dGrid< elType, Comm >
  ::mark ( int ref, const typename Traits::template Codim< 0 >::Entity &entity )
  {
    bool marked = (this->getRealImplementation( entity )).mark(ref);
    if(marked)
    {
      if(ref > 0) ++refineMarked_;
      if(ref < 0) ++coarsenMarked_;
    }
    return marked;
  }


  // get Mark of given entity
  template< ALU3dGridElementType elType, class Comm >
  inline int ALU3dGrid< elType, Comm >
  ::getMark ( const typename Traits::template Codim< 0 >::Entity &entity ) const
  {
    return this->getRealImplementation( entity ).getMark();
  }


  // global refine
  template< ALU3dGridElementType elType, class Comm >
  template< class GridImp, class DataHandle >
  inline
  void ALU3dGrid< elType, Comm >
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


  // adapt grid
  // --adapt
  template< ALU3dGridElementType elType, class Comm >
  template< class GridImp, class DataHandle >
  inline
  bool ALU3dGrid< elType, Comm >
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



  //*****************************************************************

  template< ALU3dGridElementType elType, class Comm >
  struct ALU3dGridCommHelper;

  template< ALU3dGridElementType elType >
  struct ALU3dGridCommHelper< elType, No_Comm >
  {
    typedef ALU3dGrid< elType, No_Comm > Grid;

    static bool loadBalance ( Grid &grid ) { return false; }

    template< class DataHandle >
    static bool loadBalance ( Grid &grid, DataHandle &data ) { return false; }

    template< class DataHandle, class DataType >
    static void communicate ( const Grid &grid,
                              const CommDataHandleIF< DataHandle, DataType > &data,
                              const InterfaceType iftype,
                              const CommunicationDirection dir,
                              const int level )
    {}

    template< class DataHandle, class DataType >
    static void communicate ( const Grid &grid,
                              const CommDataHandleIF< DataHandle, DataType > &data,
                              const InterfaceType iftype,
                              const CommunicationDirection dir )
    {}
  }; // ALU3dGridCommHelper

#if ALU3DGRID_PARALLEL
  template< ALU3dGridElementType elType >
  struct ALU3dGridCommHelper< elType, MPI_Comm >
  {
    typedef ALU3dGrid< elType, MPI_Comm > Grid;
    typedef ALU3DSPACE GatherScatter GatherScatterType;

    static bool loadBalance ( Grid &grid )
    {
      if( grid.comm().size() <= 1 )
        return false;

      const bool changed = grid.myGrid().duneLoadBalance();
      if( changed )
      {
        // some exchanges on ALUGrid side
        grid.myGrid().duneExchangeDynamicState();

        // calculate new maxlevel
        // reset size and things
        grid.updateStatus();

        // build new Id Set. Only do that after updateStatus, because here
        // the item lists are needed
        if( grid.globalIdSet_ )
          grid.globalIdSet_->updateIdSet();

        // unset all leaf markers
        grid.postAdapt();
      }

      return changed;
    }


    template< class DataHandle >
    static bool loadBalance ( Grid &grid, DataHandle &data )
    {
      if( grid.comm().size() <= 1 )
        return false;

      typedef typename Grid :: EntityObject EntityObject;
      typedef typename EntityObject::ImplementationType EntityImp;
      EntityObject en     ( EntityImp( grid, grid.maxLevel()) );
      EntityObject father ( EntityImp( grid, grid.maxLevel()) );
      EntityObject son    ( EntityImp( grid, grid.maxLevel()) );

      typedef ALU3DSPACE LoadBalanceElementCount< Grid, DataHandle > LDBElCountType;

      // elCount is the adaption restPro operator used during the refinement
      // cause be creating new elements on processors
      LDBElCountType elCount( grid,
                              father, Grid::getRealImplementation( father ),
                              son, Grid::getRealImplementation( son ),
                              data );

      ALU3DSPACE GatherScatterLoadBalance< Grid, DataHandle, LDBElCountType >
      gs( grid, en, Grid::getRealImplementation( en ), data, elCount );

      // call load Balance
      const bool changed = grid.myGrid().duneLoadBalance( gs, elCount );

      if( changed )
      {
        // exchange some data for internal useage
        grid.myGrid().duneExchangeDynamicState();

        // calculate new maxlevel
        // reset size and things
        grid.updateStatus();

        // build new Id Set. Only do that after updateStatus, because here
        // the item lists are needed
        if( grid.globalIdSet_ )
          grid.globalIdSet_->updateIdSet();

        // compress data, wrapper for dof manager
        gs.compress();

        grid.postAdapt();
      }
      return changed;
    }


    template< class DataHandle, class DataType >
    static void communicate ( const Grid &grid,
                              CommDataHandleIF< DataHandle, DataType > &data,
                              const InterfaceType iftype,
                              const CommunicationDirection dir,
                              const int level )
    {
      typedef CommDataHandleIF< DataHandle, DataType > DataHandleType;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 3 >::Entity > VertexObject;
      typedef typename VertexObject::ImplementationType VertexImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 2 >::Entity > EdgeObject;
      typedef typename EdgeObject::ImplementationType EdgeImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 1 >::Entity > FaceObject;
      typedef typename FaceObject::ImplementationType FaceImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 0 >::Entity> ElementObject;
      typedef typename ElementObject::ImplementationType ElementImp;

      if( grid.comm().size() > 1 )
      {
        // for level communication the level index set is needed.
        // if non-existent, then create for communicaton
        const typename Grid::LevelIndexSetImp *levelISet;
        if( !grid.levelIndexVec_[ level ] )
          levelISet = new typename Grid::LevelIndexSetImp(
            grid,
            grid.template lbegin<0>( level ),
            grid.template lend<0>( level ), level );
        else
          levelISet = grid.levelIndexVec_[ level ];

        VertexObject vx( VertexImp( grid, level ) );
        ALU3DSPACE GatherScatterLevelData< Grid, DataHandleType, 3 >
        vertexData( grid, vx, Grid::getRealImplementation( vx ), data, *levelISet, level );

        EdgeObject edge( EdgeImp( grid, level ) );
        ALU3DSPACE GatherScatterLevelData< Grid, DataHandleType, 2 >
        edgeData( grid, edge, Grid::getRealImplementation( edge ), data, *levelISet, level );

        FaceObject face( FaceImp( grid, level ) );
        ALU3DSPACE GatherScatterLevelData< Grid, DataHandleType, 1 >
        faceData( grid, face, Grid::getRealImplementation( face ), data, *levelISet, level );

        ElementObject element( ElementImp( grid, level ) );
        ALU3DSPACE GatherScatterLevelData< Grid, DataHandleType, 0 >
        elementData( grid, element, Grid::getRealImplementation( element ), data, *levelISet, level );

        doCommunication( grid, vertexData, edgeData, faceData, elementData, iftype, dir );

        if( !grid.levelIndexVec_[ level ] )
          delete levelISet;
      }
    }

    template< class DataHandle, class DataType >
    static void communicate ( const Grid &grid,
                              CommDataHandleIF< DataHandle, DataType > &data,
                              const InterfaceType iftype,
                              const CommunicationDirection dir )
    {
      typedef CommDataHandleIF< DataHandle, DataType > DataHandleType;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 3 >::Entity > VertexObject;
      typedef typename VertexObject::ImplementationType VertexImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 2 >::Entity > EdgeObject;
      typedef typename EdgeObject::ImplementationType EdgeImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 1 >::Entity > FaceObject;
      typedef typename FaceObject::ImplementationType FaceImp;
      typedef MakeableInterfaceObject< typename Grid::Traits::template Codim< 0 >::Entity> ElementObject;
      typedef typename ElementObject::ImplementationType ElementImp;

      if( grid.comm().size() > 1 )
      {
        VertexObject vx( VertexImp( grid, grid.maxLevel() ) );
        ALU3DSPACE GatherScatterLeafData< Grid, DataHandleType, 3 >
        vertexData( grid, vx, Grid::getRealImplementation( vx ), data );

        EdgeObject edge( EdgeImp( grid, grid.maxLevel() ) );
        ALU3DSPACE GatherScatterLeafData< Grid, DataHandleType, 2 >
        edgeData( grid, edge, Grid::getRealImplementation( edge ), data );

        FaceObject face( FaceImp( grid, grid.maxLevel()) );
        ALU3DSPACE GatherScatterLeafData< Grid, DataHandleType, 1 >
        faceData( grid, face, Grid::getRealImplementation( face ), data );

        ElementObject element( ElementImp( grid, grid.maxLevel() ) );
        ALU3DSPACE GatherScatterLeafData< Grid, DataHandleType, 0 >
        elementData( grid, element, Grid::getRealImplementation( element ), data );

        doCommunication( grid, vertexData, edgeData, faceData, elementData, iftype, dir );
      }
    }

    static void
    doCommunication ( const Grid &grid,
                      GatherScatterType &vertexData, GatherScatterType &edgeData,
                      GatherScatterType &faceData, GatherScatterType &elementData,
                      InterfaceType iftype, CommunicationDirection dir )
    {
      // check interface types
      if( (iftype == Overlap_OverlapFront_Interface) || (iftype == Overlap_All_Interface) )
      {
        dverb << "ALUGrid contains no overlap, therefore no communication for" << std::endl;
        dverb << "Overlap_OverlapFront_Interface or Overlap_All_Interface interfaces!" << std::endl;
      }
      // communication from border to border
      else if( iftype == InteriorBorder_InteriorBorder_Interface )
        grid.myGrid().borderBorderCommunication(vertexData,edgeData,faceData,elementData);
      // communication from interior to ghost including border
      else if( iftype == InteriorBorder_All_Interface )
      {
        if( dir == ForwardCommunication )
          grid.myGrid().interiorGhostCommunication(vertexData,edgeData,faceData,elementData);
        // reverse communiction interface (here All_InteriorBorder)
        else if( dir == BackwardCommunication )
          grid.myGrid().ghostInteriorCommunication(vertexData,edgeData,faceData,elementData);
      }
      // communication from interior to ghost including border
      else if( iftype == All_All_Interface )
        grid.myGrid().allAllCommunication(vertexData,edgeData,faceData,elementData);
      else
        DUNE_THROW( GridError, "Wrong set of parameters in ALUGridCommHelper::doCommunication" );
    }
  }; // ALU3dGridCommHelper
#endif // #if ALU3DGRID_PARALLEL



  template< ALU3dGridElementType elType, class Comm >
  inline bool ALU3dGrid< elType, Comm >::loadBalance ()
  {
    return ALU3dGridCommHelper< elType, Comm >::loadBalance( *this );
  }


  // load balance grid
  template< ALU3dGridElementType elType, class Comm >
  template< class DataHandle >
  inline bool ALU3dGrid< elType, Comm >::loadBalance ( DataHandle &data )
  {
    return ALU3dGridCommHelper< elType, Comm >::loadBalance( *this, data );
  }


  // communicate level data
  template< ALU3dGridElementType elType, class Comm >
  template <class DataHandleImp,class DataType>
  inline void ALU3dGrid< elType, Comm >::
  communicate (CommDataHandleIF<DataHandleImp,DataType> &data,
               InterfaceType iftype, CommunicationDirection dir, int level ) const
  {
    ALU3dGridCommHelper< elType, Comm >::communicate( *this, data, iftype, dir, level );
  }


  // communicate data
  template< ALU3dGridElementType elType, class Comm >
  template <class DataHandleImp, class DataType>
  inline void ALU3dGrid< elType, Comm >::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir) const
  {
    ALU3dGridCommHelper< elType, Comm >::communicate( *this, data, iftype, dir );
  }


  // return Grid name
  template< ALU3dGridElementType elType, class Comm >
  inline std::string ALU3dGrid< elType, Comm >::name ()
  {
    if( elType == hexa )
      return "ALUCubeGrid";
    else
      return "ALUSimplexGrid";
  }

} // end namespace Dune
