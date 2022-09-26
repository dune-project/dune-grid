// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include "geometry.cc"
#include "entity.cc"
#include "intersection.cc"

namespace Dune
{

  namespace Alberta
  {
    static void *adaptationDataHandler_;
  }


  template< int dim, int dimworld >
  static void checkAlbertaDimensions ()
  {
    // If this check fails, define ALBERTA_DIM accordingly
    static_assert((dimworld == Alberta::dimWorld),
                  "Template Parameter dimworld does not match "
                  "ALBERTA's DIM_OF_WORLD setting.");
  }


  // AlbertaGrid
  // -----------

  template< int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::AlbertaGrid ()
    : mesh_(),
      maxlevel_( 0 ),
      numBoundarySegments_( 0 ),
      hIndexSet_( dofNumbering_ ),
      idSet_( hIndexSet_ ),
      levelIndexVec_( (size_t)MAXL, 0 ),
      leafIndexSet_( 0 ),
      sizeCache_( *this ),
      leafMarkerVector_( dofNumbering_ ),
      levelMarkerVector_( (size_t)MAXL, MarkerVector( dofNumbering_ ) )
  {
    checkAlbertaDimensions< dim, dimworld>();
  }


  template< int dim, int dimworld >
  template< class Proj, class Impl >
  inline AlbertaGrid< dim, dimworld >
  ::AlbertaGrid ( const Alberta::MacroData< dimension> &macroData,
                  const Alberta::ProjectionFactoryInterface< Proj, Impl > &projectionFactory )
    : mesh_(),
      maxlevel_( 0 ),
      numBoundarySegments_( 0 ),
      hIndexSet_( dofNumbering_ ),
      idSet_( hIndexSet_ ),
      levelIndexVec_( (size_t)MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_( *this ),
      leafMarkerVector_( dofNumbering_ ),
      levelMarkerVector_( (size_t)MAXL, MarkerVector( dofNumbering_ ) )
  {
    checkAlbertaDimensions< dim, dimworld >();

    numBoundarySegments_ = mesh_.create( macroData, projectionFactory );
    if( !mesh_ )
      DUNE_THROW( AlbertaError, "Invalid macro data structure." );

    setup();
    hIndexSet_.create();

    calcExtras();
  }


  template< int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >
  ::AlbertaGrid ( const Alberta::MacroData< dimension> &macroData,
                  const std::shared_ptr< DuneBoundaryProjection< dimensionworld > > &projection )
    : mesh_(),
      maxlevel_( 0 ),
      numBoundarySegments_( 0 ),
      hIndexSet_( dofNumbering_ ),
      idSet_( hIndexSet_ ),
      levelIndexVec_( (size_t)MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_( *this ),
      leafMarkerVector_( dofNumbering_ ),
      levelMarkerVector_( (size_t)MAXL, MarkerVector( dofNumbering_ ) )
  {
    checkAlbertaDimensions< dim, dimworld >();

    if( projection != 0 )
    {
      Alberta::DuneGlobalBoundaryProjectionFactory< dimension > projectionFactory( projection );
      numBoundarySegments_ = mesh_.create( macroData, projectionFactory );
    }
    else
      numBoundarySegments_ = mesh_.create( macroData );
    if( !mesh_ )
      DUNE_THROW( AlbertaError, "Invalid macro data structure." );

    setup();
    hIndexSet_.create();

    calcExtras();
  }


  template < int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >
  ::AlbertaGrid ( const std::string &macroGridFileName )
    : mesh_(),
      maxlevel_( 0 ),
      hIndexSet_( dofNumbering_ ),
      idSet_( hIndexSet_ ),
      levelIndexVec_( (size_t)MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      sizeCache_( *this ),
      leafMarkerVector_( dofNumbering_ ),
      levelMarkerVector_( (size_t)MAXL, MarkerVector( dofNumbering_ ) )
  {
    checkAlbertaDimensions< dim, dimworld >();

    numBoundarySegments_ = mesh_.create( macroGridFileName );
    if( !mesh_ )
    {
      DUNE_THROW( AlbertaIOError,
                  "Grid file '" << macroGridFileName
                                << "' is not in ALBERTA macro triangulation format." );
    }

    setup();
    hIndexSet_.create();

    calcExtras();

    std::cout << typeName() << " created from macro grid file '"
              << macroGridFileName << "'." << std::endl;
  }


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::setup ()
  {
    dofNumbering_.create( mesh_ );

    levelProvider_.create( dofNumbering_ );

#if DUNE_ALBERTA_CACHE_COORDINATES
    coordCache_.create( dofNumbering_ );
#endif
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
    hIndexSet_.release();
    levelProvider_.release();
#if DUNE_ALBERTA_CACHE_COORDINATES
    coordCache_.release();
#endif
    dofNumbering_.release();

    sizeCache_.reset();

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
    typedef AlbertaGridLevelIterator< codim, pitype, const This > LevelIteratorImp;
    assert( level >= 0 );

    if( level > maxlevel_ )
      return lend< codim, pitype >( level );
    MarkerVector &markerVector = levelMarkerVector_[ level ];

    if( (codim > 0) && !markerVector.up2Date() )
      markerVector.template markSubEntities< 1 >( lbegin< 0 >( level ), lend< 0 >( level ) );

    return LevelIteratorImp( *this, &markerVector, level );
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition< pitype >::LevelIterator
  AlbertaGrid < dim, dimworld >::lend ( int level ) const
  {
    typedef AlbertaGridLevelIterator< codim, pitype, const This > LevelIteratorImp;
    assert( level >= 0 );

    return LevelIteratorImp( *this, level );
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
    typedef AlbertaGridLeafIterator< codim, pitype, const This > LeafIteratorImp;

    MarkerVector &markerVector = leafMarkerVector_;
    const int firstMarkedCodim = 2;
    if( (codim >= firstMarkedCodim) && !markerVector.up2Date() )
      markerVector.template markSubEntities< firstMarkedCodim >( leafbegin< 0 >(), leafend< 0 >() );

    return LeafIteratorImp( *this, &markerVector, maxlevel_ );
  }


  template< int dim, int dimworld >
  template< int codim, PartitionIteratorType pitype >
  inline typename AlbertaGrid< dim, dimworld >::Traits
  ::template Codim< codim >::template Partition< pitype >::LeafIterator
  AlbertaGrid< dim, dimworld >::leafend () const
  {
    typedef AlbertaGridLeafIterator< codim, pitype, const This > LeafIteratorImp;
    return LeafIteratorImp( *this, maxlevel_ );
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


  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld >::globalRefine ( int refCount )
  {
    typedef typename Traits::template Codim< 0 >::LeafIterator LeafIterator;

    // only MAXL levels allowed
    assert( (refCount >= 0) && (refCount + maxlevel_ < MAXL) );

    for( int i = 0; i < refCount; ++i )
    {
      // mark all interior elements
      const LeafIterator endit = leafend< 0 >();
      for( LeafIterator it = leafbegin< 0 >(); it != endit; ++it )
        mark( 1, *it );

      preAdapt();
      adapt();
      postAdapt();
    }
  }


  template< int dim, int dimworld >
  template< class DataHandle >
  inline void AlbertaGrid< dim, dimworld >
  ::globalRefine ( int refCount, AdaptDataHandleInterface< This, DataHandle > &handle )
  {
    typedef typename Traits::template Codim< 0 >::LeafIterator LeafIterator;

    // only MAXL levels allowed
    assert( (refCount >= 0) && (refCount + maxlevel_ < MAXL) );

    for( int i = 0; i < refCount; ++i )
    {
      // mark all interior elements
      const LeafIterator endit = leafend< 0 >();
      for( LeafIterator it = leafbegin< 0 >(); it != endit; ++it )
        mark( 1, *it );

      adapt( handle );
    }
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >::preAdapt ()
  {
    adaptationState_.preAdapt();
    return adaptationState_.coarsen();
  }


  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::postAdapt ()
  {
#ifndef NDEBUG
    if( leafIndexSet_ != 0 )
    {
      bool consistent = true;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        if( int(leafIndexSet_->size( codim )) == mesh_.size( codim ) )
          continue;
        std::cerr << "Incorrect size of leaf index set for codimension "
                  << codim << "!" << std::endl;
        std::cerr << "DUNE index set reports: " << leafIndexSet_->size( codim )
                  << std::endl;
        std::cerr << "ALBERTA mesh reports: " << mesh_.size( codim ) << std::endl;
        consistent = false;
      }
      if( !consistent )
        abort();
    }
#endif

    levelProvider_.markAllOld();
    adaptationState_.postAdapt();
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::mark( int refCount, const typename Traits::template Codim< 0 >::Entity &e )
  {
    // if not leaf entity, leave method
    if( !e.isLeaf() )
      return false;

    // don't mark macro elements for coarsening
    if( refCount < -e.level() )
      return false;

    // take back previous marking
    adaptationState_.unmark( getMark( e ) );

    // set new marking
    adaptationState_.mark( refCount );
    e.impl().elementInfo().setMark( refCount );

    return true;
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >
  ::getMark( const typename Traits::template Codim< 0 >::Entity &e ) const
  {
    return e.impl().elementInfo().getMark();
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >::adapt ()
  {
    // this is already done in postAdapt
    //levelProvider_.markAllOld();

    // adapt mesh
    hIndexSet_.preAdapt();
    const bool refined = mesh_.refine();
    const bool coarsened = (adaptationState_.coarsen() ? mesh_.coarsen() : false);
    adaptationState_.adapt();
    hIndexSet_.postAdapt();

    if( refined || coarsened )
      calcExtras();

    // return true if elements were created
    return refined;
  }


  template< int dim, int dimworld >
  template< class DataHandle >
  inline bool AlbertaGrid < dim, dimworld >
  ::adapt ( AdaptDataHandleInterface< This, DataHandle > &handle )
  {
    preAdapt();

    typedef Alberta::AdaptRestrictProlongHandler
    < This, AdaptDataHandleInterface< This, DataHandle > >
    DataHandler;
    DataHandler dataHandler( *this, handle );

    typedef AdaptationCallback< DataHandler > Callback;
    typename Callback::DofVectorPointer callbackVector;
    callbackVector.create( dofNumbering_.emptyDofSpace(), "Adaptation Callback" );
    callbackVector.template setupInterpolation< Callback >();
    callbackVector.template setupRestriction< Callback >();
    if( Callback::DofVectorPointer::supportsAdaptationData )
      callbackVector.setAdaptationData( &dataHandler );
    else
      Alberta::adaptationDataHandler_ = &dataHandler;

    bool refined = adapt();

    if( !Callback::DofVectorPointer::supportsAdaptationData )
      Alberta::adaptationDataHandler_ = 0;
    callbackVector.release();

    postAdapt();
    return refined;
  }


  template< int dim, int dimworld >
  inline const Alberta::GlobalVector &
  AlbertaGrid< dim, dimworld >
  ::getCoord ( const ElementInfo &elementInfo, int vertex ) const
  {
    assert( (vertex >= 0) && (vertex <= dim) );
#if DUNE_ALBERTA_CACHE_COORDINATES
    return coordCache_( elementInfo, vertex );
#else
    return elementInfo.coordinate( vertex );
#endif
  }


  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::maxLevel() const
  {
    return maxlevel_;
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >::size ( int level, int codim ) const
  {
    return ((level >= 0) && (level <= maxlevel_) ? sizeCache_.size( level, codim ) : 0);
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >::size ( int level, GeometryType type ) const
  {
    return ((level >= 0) && (level <= maxlevel_) ? sizeCache_.size( level, type ) : 0);
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >::size ( int codim ) const
  {
    assert( sizeCache_.size( codim ) == mesh_.size( codim ) );
    return mesh_.size( codim );
  }


  template< int dim, int dimworld >
  inline int AlbertaGrid< dim, dimworld >::size ( GeometryType type ) const
  {
    return sizeCache_.size( type );
  }


  template < int dim, int dimworld >
  inline const typename AlbertaGrid < dim, dimworld > :: Traits :: LevelIndexSet &
  AlbertaGrid < dim, dimworld > :: levelIndexSet (int level) const
  {
    // assert that given level is in range
    assert( (level >= 0) && (level < (int)levelIndexVec_.size()) );

    if( levelIndexVec_[ level ] == 0 )
    {
      levelIndexVec_[ level ] = new typename GridFamily::LevelIndexSetImp ( dofNumbering_ );
      levelIndexVec_[ level ]->update( lbegin< 0 >( level ), lend< 0 >( level ) );
    }
    return *(levelIndexVec_[ level ]);
  }

  template < int dim, int dimworld >
  inline const typename AlbertaGrid < dim, dimworld > :: Traits :: LeafIndexSet &
  AlbertaGrid < dim, dimworld > :: leafIndexSet () const
  {
    if( leafIndexSet_ == 0 )
    {
      leafIndexSet_ = new typename GridFamily::LeafIndexSetImp( dofNumbering_ );
      leafIndexSet_->update( leafbegin< 0 >(), leafend< 0 >() );
    }
    return *leafIndexSet_;
  }


  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::calcExtras ()
  {
    // determine new maxlevel
    maxlevel_ = levelProvider_.maxLevel();
    assert( (maxlevel_ >= 0) && (maxlevel_ < MAXL) );

    // unset up2Dat status, if lbegin is called then this status is updated
    for( int l = 0; l < MAXL; ++l )
      levelMarkerVector_[ l ].clear();

    // unset up2Dat status, if leafbegin is called then this status is updated
    leafMarkerVector_.clear();

    sizeCache_.reset();

    // update index sets (if they exist)
    if( leafIndexSet_ != 0 )
      leafIndexSet_->update( leafbegin< 0 >(), leafend< 0 >() );
    for( unsigned int level = 0; level < levelIndexVec_.size(); ++level )
    {
      if( levelIndexVec_[ level ] )
        levelIndexVec_[ level ]->update( lbegin< 0 >( level ), lend< 0 >( level ) );
    }
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::writeGrid ( const std::string &filename, ctype time ) const
  {
    if( filename.size() <= 0 )
      DUNE_THROW( AlbertaIOError, "No filename given to writeGridXdr." );
    return (mesh_.write( filename, time ) && hIndexSet_.write( filename ));
  }


  template< int dim, int dimworld >
  inline bool AlbertaGrid< dim, dimworld >
  ::readGrid ( const std::string &filename, ctype &time )
  {
    //removeMesh();

    if( filename.size() <= 0 )
      return false;

    numBoundarySegments_ = mesh_.read( filename, time );
    if( !mesh_ )
      DUNE_THROW( AlbertaIOError, "Could not read grid file: " << filename << "." );

    setup();
    hIndexSet_.read( filename );

    // calc maxlevel and indexOnLevel and so on
    calcExtras();

    return true;
  }


  // AlbertaGrid::AdaptationCallback
  // -------------------------------

  template< int dim, int dimworld >
  template< class DataHandler >
  struct AlbertaGrid< dim, dimworld >::AdaptationCallback
  {
    static const int dimension = dim;

    typedef Alberta::DofVectorPointer< Alberta::GlobalVector > DofVectorPointer;
    typedef Alberta::Patch< dimension > Patch;

    static DataHandler &getDataHandler ( const DofVectorPointer &dofVector )
    {
      DataHandler *dataHandler;
      if( DofVectorPointer::supportsAdaptationData )
        dataHandler = dofVector.getAdaptationData< DataHandler >();
      else
        dataHandler = (DataHandler *)Alberta::adaptationDataHandler_;
      assert( dataHandler != 0 );
      return *dataHandler;
    }

    static void interpolateVector ( const DofVectorPointer &dofVector,
                                    const Patch &patch )
    {
      DataHandler &dataHandler = getDataHandler( dofVector );
      for( int i = 0; i < patch.count(); ++i )
        dataHandler.prolongLocal( patch, i );
    }

    static void restrictVector ( const DofVectorPointer &dofVector,
                                 const Patch &patch )
    {
      DataHandler &dataHandler = getDataHandler( dofVector );
      for( int i = 0; i < patch.count(); ++i )
        dataHandler.restrictLocal( patch, i );
    }
  };

} // namespace Dune

#endif
