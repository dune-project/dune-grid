// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITY_CC
#define DUNE_ALBERTA_ENTITY_CC

#include <dune/grid/albertagrid/boundary.hh>

namespace Dune
{

  //*************************************************************************
  //
  //  --AlbertaGridEntity
  //  --Entity
  //
  //*************************************************************************
  //
  //  codim > 0
  //
  // The Geometry is prescribed by the EL_INFO struct of ALBERTA MESH
  // the pointer to this struct is set and get by setElInfo and
  // getElInfo.
  //*********************************************************************8
  template<int codim, int dim, class GridImp>
  inline AlbertaGridEntity<codim,dim,GridImp>::
  AlbertaGridEntity(const GridImp &grid, int level,
                    ALBERTA TRAVERSE_STACK * travStack)
    : grid_( grid ),
      elInfo_(0),
      element_(0),
      travStack_(travStack),
      level_ ( level ),
      geo_ (GeometryImp()),
      builtgeometry_    (false),
      localFatherCoords_(),
      localFCoordCalced_(false),
      subEntity_( -1 )
  {}


  template<int codim, int dim, class GridImp>
  inline AlbertaGridEntity< codim, dim, GridImp >
  :: AlbertaGridEntity ( const This &other )
    : grid_( other.grid_ ),
      elInfo_( other.elInfo_ ),
      element_( elInfo_ != NULL ? elInfo_->el : NULL ),
      travStack_( other.travStack_ ),
      level_( other.level_ ),
      geo_( other.geo_ ),
      builtgeometry_( false ),
      localFatherCoords_(),
      localFCoordCalced_( false ),
      subEntity_( other.subEntity_ )
  {}


  template<int codim, int dim, class GridImp>
  inline void AlbertaGridEntity<codim,dim,GridImp>::
  setTraverseStack(ALBERTA TRAVERSE_STACK * travStack)
  {
    travStack_ = travStack;
  }

  template<int codim, int dim, class GridImp>
  inline AlbertaGridEntity<codim,dim,GridImp>::
  AlbertaGridEntity(const GridImp &grid, int level, bool)
    : grid_(grid),
      elInfo_(0),
      element_(0),
      travStack_(0),
      level_ (level),
      geo_ (GeometryImp()),
      builtgeometry_(false),
      localFatherCoords_(),
      localFCoordCalced_(false),
      subEntity_( -1 )
  {}

  template<int codim, int dim, class GridImp>
  inline PartitionType AlbertaGridEntity <codim,dim,GridImp>::
  partitionType () const
  {
    return InteriorEntity;
  }

  template< int codim, int dim, class GridImp >
  inline bool AlbertaGridEntity< codim, dim, GridImp >
  :: equals ( const This &other ) const
  {
    const ALBERTA EL *e2 = other.getElement();

    // if both element null then they are equal
    if( (!e2) && (!element_) )
      return true;

    return ((element_ == e2) && (subEntity_ == other.subEntity_));
  }

  template<int codim, int dim, class GridImp>
  inline ALBERTA EL_INFO* AlbertaGridEntity<codim,dim,GridImp>::
  getElInfo() const
  {
    return elInfo_;
  }

  template<int codim, int dim, class GridImp>
  inline ALBERTA EL * AlbertaGridEntity<codim,dim,GridImp>::
  getElement() const
  {
    return element_;
  }

  template<int codim, int dim, class GridImp>
  inline void AlbertaGridEntity<codim,dim,GridImp>::
  removeElInfo()
  {
    elInfo_  = 0;
    element_ = 0;
    builtgeometry_ = false;
  }

  template< int codim, int dim, class GridImp >
  inline void AlbertaGridEntity< codim, dim, GridImp >
  ::setElInfo( ALBERTA EL_INFO *elInfo, int subEntity )
  {
    elInfo_ = elInfo;
    subEntity_ = subEntity;
    element_ = (elInfo_ != NULL ? elInfo_->el : NULL );
    builtgeometry_ = geoImp().builtGeom( grid_, elInfo_, subEntity );
    localFCoordCalced_ = false;
  }

  template< int codim, int dim, class GridImp >
  inline void AlbertaGridEntity< codim, dim, GridImp >
  ::setEntity( const This &other )
  {
    setElInfo( other.elInfo_, other.subEntity_ );
    setLevel( other.level_ );
  }

  template<int codim, int dim, class GridImp>
  inline void AlbertaGridEntity<codim,dim,GridImp>::
  setLevel(int level)
  {
    level_  = level;
  }

  template<int codim, int dim, class GridImp>
  inline int AlbertaGridEntity<codim,dim,GridImp>::
  level() const
  {
    return level_;
  }

#if 0
  // default
  template< class GridImp, int codim, int cdim >
  struct AlbertaGridBoundaryId
  {
    static int boundaryId ( const ALBERTA EL_INFO *elInfo, int subEntity )
    {
      return 0;
    }
  };

  // faces in 2d and 3d
  template< class GridImp >
  struct AlbertaGridBoundaryId< GridImp, 1, 3 >
  {
    static int boundaryId ( const ALBERTA EL_INFO *elInfo, int subEntity )
    {
      return EDGE_BOUNDARY_ID( elInfo, subEntity );
    }
  };

  template< class GridImp >
  struct AlbertaGridBoundaryId< GridImp, 1, 2 >
  {
    static int boundaryId ( const ALBERTA EL_INFO *elInfo, int subEntity )
    {
      return EDGE_BOUNDARY_ID( elInfo, subEntity );
    }
  };

  // vertices in 2d and 3d
  template< class GridImp, int dim >
  struct AlbertaGridBoundaryId< GridImp, dim, dim >
  {
    static int boundaryId ( const ALBERTA EL_INFO *elInfo, int subEntity )
    {
      return VERTEX_BOUNDARY_ID( elInfo, subEntity );
    }
  };
#endif

  template< int codim, int dim, class GridImp >
  inline int AlbertaGridEntity< codim, dim, GridImp >::boundaryId () const
  {
    return Alberta::template boundaryId< dim, codim >( elInfo_, subEntity_ );
    //return AlbertaGridBoundaryId< GridImp, codim, dim >::boundaryId( elInfo_, subEntity_ );
  }

  template<int codim, int dim, class GridImp>
  inline int AlbertaGridEntity<codim,dim,GridImp>::
  getFEVnum() const
  {
    return subEntity_;
  }

  template<int cd, int dim, class GridImp>
  inline const typename AlbertaGridEntity<cd,dim,GridImp>::Geometry &
  AlbertaGridEntity<cd,dim,GridImp>::geometry() const
  {
    assert(builtgeometry_ == true);
    return geo_;
  }

  template< int cd, int dim, class GridImp >
  inline GeometryType AlbertaGridEntity< cd, dim, GridImp> :: type () const
  {
    return GeometryType( GeometryType :: simplex, mydimension );
  }

  //************************************
  //
  //  --AlbertaGridEntity codim = 0
  //  --0Entity codim = 0
  //
  template<int dim, class GridImp>
  inline AlbertaGridEntity <0,dim,GridImp>::
  AlbertaGridEntity(const GridImp &grid, int level, bool leafIt )
    : grid_(grid)
      , level_ (level)
      , travStack_ (0)
      , elInfo_ (0)
      , element_(0)
      , geoObj_( GeometryImp() )
      , geo_( grid_.getRealImplementation(geoObj_) )
      , builtgeometry_ (false)
      , leafIt_ ( leafIt )
  {}

  template<int dim, class GridImp>
  inline AlbertaGridEntity <0,dim,GridImp>::
  AlbertaGridEntity(const AlbertaGridEntity & org)
    : grid_(org.grid_)
      , level_ (org.level_)
      , travStack_ (org.travStack_)
      , elInfo_ (org.elInfo_)
      , element_( (elInfo_) ? (elInfo_->el) : 0)
      , geoObj_( org.geoObj_ )
      , geo_( grid_.getRealImplementation(geoObj_) )
      , builtgeometry_ (false)
      , leafIt_ ( org.leafIt_ )
  {}

  template<int dim, class GridImp>
  inline int AlbertaGridEntity <0,dim,GridImp>::
  boundaryId() const
  {
    // elements are always inside of our Domain
    return 0;
  }

  template<int dim, class GridImp>
  inline bool AlbertaGridEntity <0,dim,GridImp>::
  isNew () const
  {
    assert( element_ && elInfo_ );
    assert( element_ == elInfo_->el );
    return grid_.checkElNew( element_ );
  }

  template<int dim, class GridImp>
  inline bool AlbertaGridEntity <0,dim,GridImp>::
  mightVanish () const
  {
    assert( element_ && elInfo_ );
    assert( element_ == elInfo_->el );
    return ( element_->mark < 0 );
  }

  template< int dim, class GridImp >
  inline bool AlbertaGridEntity< 0, dim, GridImp >
  ::hasBoundaryIntersections () const
  {
    assert( elInfo_ != NULL );
    bool isBoundary = false;
    for( int i = 0; i < dim+1; ++i )
      isBoundary |= Alberta::isBoundary( elInfo_, i );
    return isBoundary;
  }

  template<int dim, class GridImp>
  inline PartitionType AlbertaGridEntity <0,dim,GridImp>::
  partitionType () const
  {
    return InteriorEntity;
  }

  template<int dim, class GridImp>
  inline bool AlbertaGridEntity <0,dim,GridImp>::isLeaf() const
  {
    assert( element_ && elInfo_ );
    assert( element_ == elInfo_->el );

    // if no child exists, then this element is leaf element
    return IS_LEAF_EL(element_);
  }

  //***************************

  template<int dim, class GridImp>
  inline void AlbertaGridEntity <0,dim,GridImp>::
  makeDescription()
  {
    elInfo_  = 0;
    element_ = 0;
    builtgeometry_ = false;
  }

  template<int dim, class GridImp>
  inline bool AlbertaGridEntity<0,dim,GridImp>::
  equals (const AlbertaGridEntity<0,dim,GridImp> & i) const
  {
    // compare element pointer which are unique
    return (element_ == i.getElement());
  }

  template<int dim, class GridImp>
  inline void AlbertaGridEntity <0,dim,GridImp>::
  setTraverseStack(ALBERTA TRAVERSE_STACK * travStack)
  {
    travStack_ = travStack;
  }

  //*****************************************************************
  // count
  template <class GridImp, int dim, int cc> struct AlbertaGridCount {
    static int count () { return dim+1; }
  };

  // specialisation for codim 0
  template <class GridImp, int dim> struct AlbertaGridCount<GridImp,dim,0> {
    static int count () { return 1; }
  };

  // specialisation for edges in 3d
  template <class GridImp> struct AlbertaGridCount<GridImp,3,2> {
    static int count () { return 6; }
  };

  template<int dim, class GridImp> template <int cc>
  inline int AlbertaGridEntity <0,dim,GridImp>::count () const
  {
    return AlbertaGridCount<GridImp,dim,cc>::count();
  }

  //*****************************************************************
  template< int dim, class GridImp >
  template< int codim >
  inline typename AlbertaGridEntity< 0, dim, GridImp >
  ::template Codim< codim >::EntityPointer
  AlbertaGridEntity< 0, dim, GridImp >::entity ( int i ) const
  {
    typedef AlbertaGridEntityPointer< codim, GridImp > EntityPointerImpl;
    return EntityPointerImpl( grid_, travStack_, level(), elInfo_, i );
  }

  template<int dim, class GridImp>
  inline ALBERTA EL_INFO* AlbertaGridEntity <0,dim,GridImp>::
  getElInfo() const
  {
    return elInfo_;
  }

  template<int dim, class GridImp>
  inline ALBERTA EL * AlbertaGridEntity <0,dim,GridImp>::
  getElement() const
  {
    return element_;
  }

  template<int dim, class GridImp>
  inline int AlbertaGridEntity <0,dim,GridImp>::
  level() const
  {
    return level_;
  }

  template<int dim, class GridImp>
  inline void AlbertaGridEntity<0,dim,GridImp>::
  removeElInfo()
  {
    elInfo_  = 0;
    element_ = 0;
    builtgeometry_ = false;
    level_ = -1;
  }

  template<int dim, class GridImp>
  inline void AlbertaGridEntity <0,dim,GridImp>::
  setElInfo(ALBERTA EL_INFO * elInfo, int , int , int )
  {
    // just set elInfo and element
    elInfo_ = elInfo;
    if(elInfo_)
    {
      element_ = elInfo_->el;
      level_ = grid_.getLevelOfElement( element_ );
    }
    else
    {
      level_ = -1;
      element_ = 0;
    }
    builtgeometry_ = false;
  }

  template<int dim, class GridImp>
  inline void AlbertaGridEntity<0,dim,GridImp>::
  setEntity(const AlbertaGridEntity<0,dim,GridImp> & org)
  {
    setElInfo(org.elInfo_);
    setTraverseStack(org.travStack_);
  }

  template<int dim, class GridImp>
  inline const typename AlbertaGridEntity <0,dim,GridImp>::Geometry &
  AlbertaGridEntity <0,dim,GridImp>::geometry() const
  {
    assert( elInfo_ && element_ );
    // geometry is only build on demand
    if( !builtgeometry_ )
      builtgeometry_ = geo_.builtGeom( grid_, elInfo_, 0 );
    assert( builtgeometry_ );
    return geoObj_;
  }

  template< int dim, class GridImp >
  inline GeometryType AlbertaGridEntity< 0, dim, GridImp> :: type () const
  {
    return GeometryType( GeometryType :: simplex, mydimension );
  }

  // --father
  template< int dim, class GridImp >
  inline typename AlbertaGridEntity< 0, dim, GridImp >::EntityPointer
  AlbertaGridEntity< 0, dim, GridImp >::father () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImpl;
    //std::cout << "level_ = " << level_ << "\n";
    //
    ALBERTA EL_INFO *fatherInfo
      = ALBERTA AlbertHelp::getFatherInfo(travStack_,elInfo_,level_);
    assert( fatherInfo != NULL );
    // check father element pointer
    assert( elInfo_->parent == fatherInfo->el );

    //std::cout << "Father of el[" << grid_.getElementNumber(element_) << "] is father[" << grid_.getElementNumber(fatherInfo->el) << "]\n";

    const int fatherLevel = level_ - 1;
    assert( fatherLevel >= 0 );

    //std::cout << "set father with level " << fatherLevel << " | stack used = " << travStack_->stack_used << "\n";

    assert( (fatherLevel == fatherInfo->level) );

    return EntityPointerImpl( grid_, travStack_, fatherLevel, fatherInfo, 0 );
  }

  template< int dim, class GridImp >
  inline int AlbertaGridEntity< 0, dim, GridImp > :: nChild () const
  {
    // get father and check which child we have
    const ALBERTA EL *father = elInfo_->parent;
    assert( father );

    const int child = (father->child[ 0 ] == element_ ? 0 : 1);
    assert( father->child[ child ] == element_ );
    return child;
  }

  template<>
  inline const AlbertaGridEntity <0,2,const AlbertaGrid<2,2> >::Geometry &
  AlbertaGridEntity <0,2,const AlbertaGrid<2,2> >::geometryInFather() const
  {
    typedef AlbertaGridLocalGeometryProvider< 2, 2 > LocalGeoProvider;
    return LocalGeoProvider::instance().geometryInFather( nChild() );
  }

  template<>
  inline const AlbertaGridEntity <0,3,const AlbertaGrid<3,3> >::Geometry &
  AlbertaGridEntity <0,3,const AlbertaGrid<3,3> >::geometryInFather() const
  {
    // see Alberta Docu for definition  of el_type, values are 0,1,2
    const int orientation =
#if DIM == 3
      (elInfo_->el_type == 1) ? -1 :
#endif
      1;
    typedef AlbertaGridLocalGeometryProvider< 3, 3 > LocalGeoProvider;
    return LocalGeoProvider::instance().geometryInFather( nChild(), orientation );
  }

} // namespace Dune

#undef ALBERTA_CHAR
#endif
