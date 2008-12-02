// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_CC
#define DUNE_ALBERTA_INTERSECTION_CC

#include <dune/grid/albertagrid/intersection.hh>

namespace Dune
{

  template< class GridImp >
  inline AlbertaGridIntersectionIterator< GridImp >
  ::AlbertaGridIntersectionIterator ( const GridImp &grid, int level )
    : grid_(grid),
      level_ (level),
      neighborCount_( dimension+1 ),
      elInfo_ (0),
      fakeNeighObj_(LocalGeometryImp()),
      fakeSelfObj_ (LocalGeometryImp()),
      neighGlobObj_(LocalGeometryImp()),
      neighElInfo_ ()
  {}

  template< class GridImp >
  inline void
  AlbertaGridIntersectionIterator< GridImp >
  ::first ( const EntityImp &entity, int level )
  {
    level_ = level;
    neighborCount_ = 0;
    builtNeigh_    = false;
    elInfo_        = entity.getElInfo();
    this->leafIt_  = entity.leafIt();
    assert( elInfo_ );
  }

  template< class GridImp >
  inline void
  AlbertaGridIntersectionIterator< GridImp >::done ()
  {
    level_ = -1;
    neighborCount_ = dimension+1;
    builtNeigh_    = false;
    elInfo_        = NULL;
  }


  // copy constructor
  template< class GridImp >
  inline AlbertaGridIntersectionIterator< GridImp >
  ::AlbertaGridIntersectionIterator ( const This &other )
    : grid_( other.grid_ ),
      level_( other.level_ ),
      neighborCount_( other.neighborCount_ ),
      builtNeigh_( false ),
      leafIt_( other.leafIt_ ),
      elInfo_ ( other.elInfo_ ),
      fakeNeighObj_( LocalGeometryImp() ),
      fakeSelfObj_ ( LocalGeometryImp() ),
      neighGlobObj_( LocalGeometryImp() ),
      neighElInfo_()
  {}


  // assignment operator
  template< class GridImp >
  inline void
  AlbertaGridIntersectionIterator< GridImp >::assign ( const This &other )
  {
    // only assign iterators from the same grid
    assert( &(this->grid_) == &(other.grid_) );

    level_ = other.level_;
    neighborCount_ = other.neighborCount_;
    elInfo_ = other.elInfo_;
    builtNeigh_ = false;
    leafIt_ = other.leafIt_;
  }


  template< class GridImp >
  inline bool
  AlbertaGridIntersectionIterator< GridImp >::equals ( const This &other ) const
  {
    const ALBERTA EL *e1 = (elInfo_ != NULL ? elInfo_->el : NULL);
    const ALBERTA EL *e2 = (other.elInfo_ != NULL ? other.elInfo_->el : NULL);
    return (e1 == e2) && (neighborCount_ == other.neighborCount_);
  }

  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>::increment()
  {
    builtNeigh_ = false;
    // is like go to the next neighbour
    ++neighborCount_;

    // (dim+1) is neigbourCount for end iterators
    if( neighborCount_ > dimension )
    {
      this->done();
      return;
    }

    /*
       // check whether neighbor has same level
       if( (this->neighbor()) && (!this->leafIt()) )
       {
       // if levels are not equal go next
       if( !neighborHasSameLevel () ) increment ();
       }
     */
  }

  template< class GridImp >
  inline typename AlbertaGridIntersectionIterator< GridImp >::EntityPointer
  AlbertaGridIntersectionIterator< GridImp >::outside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImp;

    if( !builtNeigh_ )
    {
      assert( elInfo_ );
      // just copy elInfo and then set some values
      std::memcpy( &neighElInfo_, elInfo_, sizeof( ALBERTA EL_INFO ) );

      setupVirtEn();

      assert( level_ == elInfo_->level );
      assert( leafIt() || (elInfo_->level == neighElInfo_.level) );
    }

    assert( builtNeigh_ );
    assert( neighElInfo_.el );
    const int neighborLevel = grid_.getLevelOfElement( neighElInfo_.el );
    return EntityPointerImp( grid_, neighborLevel, &neighElInfo_, 0 );
  }

  template< class GridImp >
  inline typename AlbertaGridIntersectionIterator< GridImp >::EntityPointer
  AlbertaGridIntersectionIterator< GridImp >::inside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImp;
    assert( elInfo_ != NULL );
    return EntityPointerImp( grid_, (int)elInfo_->level, elInfo_, 0 );
  }

  template< class GridImp >
  inline int
  AlbertaGridIntersectionIterator<GridImp>::boundaryId () const
  {
    // id of interior intersections is 0
    if( ! boundary() ) return 0;
    assert(elInfo_);
    assert( WALL_BOUND(elInfo_, dimension, neighborCount_ ) != 0 );
    return WALL_BOUND(elInfo_, dimension, neighborCount_ );
  }

  template< class GridImp >
  inline bool AlbertaGridIntersectionIterator<GridImp>::conforming() const
  {
    return true;
  }

  template< class GridImp >
  inline bool AlbertaGridIntersectionIterator< GridImp >::boundary() const
  {
    assert( elInfo_ != NULL );
    return Alberta::isBoundary( elInfo_, neighborCount_ );
  }

  template< class GridImp >
  inline bool AlbertaGridIntersectionIterator<GridImp>::neighbor() const
  {
    // use ALBERTA macro to get neighbour
    assert( elInfo_ );
    return (NEIGH(elInfo_->el,elInfo_)[neighborCount_] != 0);
  }

  template<class GridImp>
  inline const typename AlbertaGridIntersectionIterator< GridImp >::NormalVector
  AlbertaGridIntersectionIterator< GridImp >
  ::unitOuterNormal ( const LocalCoordType &local ) const
  {
    NormalVector normal;
    calcOuterNormal( normal );
    normal *= (1.0 / normal.two_norm());
    return normal;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator< GridImp >::NormalVector
  AlbertaGridIntersectionIterator< GridImp >
  ::integrationOuterNormal ( const LocalCoordType &local ) const
  {
    NormalVector normal;
    calcOuterNormal( normal );
    return normal;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator< GridImp >::NormalVector
  AlbertaGridIntersectionIterator< GridImp >
  ::outerNormal( const LocalCoordType &local ) const
  {
    return integrationOuterNormal( local );
  }

  template< class GridImp >
  inline void
  AlbertaGridIntersectionIterator<GridImp>::calcOuterNormal( NormalVector &n ) const
  {
    DUNE_THROW( NotImplemented, "AlbertaGrid: outer normal for dim != dimworld "
                "has not been implemented, yet." );
  }


  template<>
  inline void
  AlbertaGridIntersectionIterator< const AlbertaGrid< 2, 2 > >
  ::calcOuterNormal ( NormalVector &n ) const
  {
    assert( elInfo_ );
    const ALBERTA REAL_D & coordOne = grid_.getCoord(elInfo_,(neighborCount_+1)%3);
    const ALBERTA REAL_D & coordTwo = grid_.getCoord(elInfo_,(neighborCount_+2)%3);

    n[ 0 ] = -(coordOne[ 1 ] - coordTwo[ 1 ]);
    n[ 1 ] =   coordOne[ 0 ] - coordTwo[ 0 ];
  }


  template<>
  inline void
  AlbertaGridIntersectionIterator< const AlbertaGrid< 3, 3 > >
  ::calcOuterNormal ( NormalVector &n ) const
  {
    assert( elInfo_ != NULL );

    // in this case the orientation is negative, multiply by -1
#if (DUNE_ALBERTA_VERSION >= 0x200) || (DIM == 3)
    const albertCtype val = (elInfo_->orientation > 0) ? 1.0 : -1.0;
#else
    const albertCtype val = 1.0;
#endif

    // neighborCount_ is the local face number
    const int *localFaces = ALBERTA AlbertHelp::localAlbertaFaceNumber[ neighborCount_ ];

    const ALBERTA REAL_D &coordZero = grid_.getCoord( elInfo_, localFaces[0] );
    const ALBERTA REAL_D &coordOne  = grid_.getCoord( elInfo_, localFaces[1] );
    const ALBERTA REAL_D &coordTwo  = grid_.getCoord( elInfo_, localFaces[2] );

    FieldVector< ctype, dimensionworld > u;
    FieldVector< ctype, dimensionworld > v;
    for( int i = 0; i < dimension; ++i )
    {
      v[ i ] = coordOne[ i ] - coordZero[ i ];
      u[ i ] = coordTwo[ i ] - coordOne[ i ];
    }

    // outNormal_ has length 3
    for( int i = 0; i < dimension; ++i )
    {
      const int j = (i+1)%dimension;
      const int k = (i+2)%dimension;
      n[ i ] = val * (u[ j ] * v[ k ] - u[ k ] * v[ j ]);
    }
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator< GridImp >::LocalGeometry &
  AlbertaGridIntersectionIterator< GridImp >::intersectionSelfLocal () const
  {
    LocalGeometryImp &geo = GridImp::getRealImplementation( fakeSelfObj_ );
    if( !geo.builtLocalGeom( inside()->geometry(), intersectionGlobal(), elInfo_, neighborCount_ ) )
      DUNE_THROW( AlbertaError, "internal error in intersectionSelfLocal." );
    return fakeSelfObj_;
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator< GridImp >::LocalGeometry &
  AlbertaGridIntersectionIterator< GridImp >::intersectionNeighborLocal () const
  {
    assert( neighbor() );
    LocalGeometryImp &geo = GridImp::getRealImplementation( fakeNeighObj_ );
    if( !geo.builtLocalGeom( outside()->geometry(), intersectionGlobal(), &neighElInfo_, neighborCount_ ) )
      DUNE_THROW( AlbertaError, "internal error in intersectionNeighborLocal." );
    return fakeNeighObj_;
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator< GridImp >::Geometry &
  AlbertaGridIntersectionIterator< GridImp >::intersectionGlobal () const
  {
    assert( elInfo_ != NULL );

    GeometryImp &geo = GridImp::getRealImplementation( neighGlobObj_ );
    if( !geo.builtGeom( grid_, elInfo_, neighborCount_ ) )
      DUNE_THROW( AlbertaError, "internal error in intersectionGlobal." );
    return neighGlobObj_;
  }


  template< class GridImp >
  inline int AlbertaGridIntersectionIterator<GridImp>::
  level () const
  {
    assert( level_ >= 0 );
    return level_;
  }

  template< class GridImp >
  inline int AlbertaGridIntersectionIterator<GridImp>::
  numberInSelf () const
  {
    return neighborCount_;
  }

  template< class GridImp >
  inline int AlbertaGridIntersectionIterator<GridImp>::
  numberInNeighbor () const
  {
    assert( elInfo_ );
    return elInfo_->opp_vertex[neighborCount_];
  }

  template <class GridImp>
  inline int AlbertaGridIntersectionIterator<GridImp>::
  twistInSelf() const
  {
    // always 0 for indside
    return 0;
  }

  template <class GridImp>
  inline int AlbertaGridIntersectionIterator<GridImp>::
  twistInNeighbor() const
  {
    return twist_;
  }

  // setup neighbor element with the information of elInfo_
  template< class GridImp >
  inline bool AlbertaGridIntersectionIterator<GridImp>::neighborHasSameLevel () const
  {
    assert( neighbor() );
    const ALBERTA EL * myEl    = elInfo_->el;
    const ALBERTA EL * neighEl = NEIGH(elInfo_->el,elInfo_)[neighborCount_];

    return this->grid_.getLevelOfElement( myEl ) ==
           this->grid_.getLevelOfElement( neighEl );
  }

  //*****************************************
  //  setup for 2d
  //*****************************************
  template <class GridImp, int dimworld , int dim >
  struct SetupVirtualNeighbour;


  #if DIM == 2
  template <class GridImp>
  struct SetupVirtualNeighbour<GridImp,2,2>
  {
    enum { dim      = 2 };
    enum { dimworld = 2 };
    static int setupNeighInfo(GridImp & grid, const ALBERTA EL_INFO * elInfo,
                              const int vx, const int nb, ALBERTA EL_INFO * neighInfo)
    {

  #ifdef CALC_COORD
      // vx is the face number in the neighbour
      const int (& neighmap)[dim] = ALBERTA AlbertHelp :: localTriangleFaceNumber[vx];
      // neighborCount is the face number in the actual element
      const int (& selfmap) [dim] = ALBERTA AlbertHelp :: localTriangleFaceNumber[nb];

      // copy the two edge vertices
      for(int i=0; i<dim; i++)
      {
        const ALBERTA REAL_D & coord = elInfo->coord[ selfmap[i] ];
        // here the twist is simple, just swap the vertices, and do this by hand
        ALBERTA REAL_D & newcoord    = neighInfo->coord[ neighmap[(dim-1) - i] ];
        for(int j=0; j<dimworld; j++) newcoord[j] = coord[j];
      }
  #endif

      //****************************************
      // twist is always 1
      return 1;
    }
  };
  #endif

  //******************************************
  //  setup for 3d
  //******************************************
  #if DIM == 3
  template <class GridImp>
  struct SetupVirtualNeighbour<GridImp,3,3>
  {
    enum { dim      = 3 };
    enum { dimworld = 3 };
    static int setupNeighInfo(GridImp & grid, const ALBERTA EL_INFO * elInfo,
                              const int vx, const int nb, ALBERTA EL_INFO * neighInfo)
    {
      // the face might be twisted when look from different elements
      // default is no, and then rthe orientation is -1
      int facemap[dim]   = {0,1,2};
      bool rightOriented = false;
      {
        int myvx[dim];
        int neighvx[dim];

        const int * vxmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[vx];
        const int * nbmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[nb];

        bool allRight = true;
        for(int i=0; i<dim; i++)
        {
          myvx[i]    = grid.getVertexNumber(elInfo->el   , nbmap[i] );
          neighvx[i] = grid.getVertexNumber(neighInfo->el, vxmap[i] );
          if( myvx[i] != neighvx[i] ) allRight = false;
        }

        // note: if the vertices are equal then the face in the neighbor
        // is not oriented right, because all face are oriented math. pos when
        // one looks from the outside of the element.
        // if the vertices are the same, the face in the neighbor is therefore
        // wrong oriented
        if( !allRight )
        {
          for(int i=0; i<dim; i++)
          {
            if(myvx[i] != neighvx[i])
            {
              for(int k=1; k<dim; k++)
              {
                int newvx = (i+k) % dim;
                if( myvx[i] == neighvx[newvx] ) facemap[i] = newvx;
              }
            }
          }
          rightOriented = true;
        }
      }

      // TODO check infulence of orientation
      // is used when calculation the outerNormal
      neighInfo->orientation = ( rightOriented ) ? elInfo->orientation : -elInfo->orientation;

  #ifdef CALC_COORD
      // vx is the face number in the neighbour
      const int * neighmap = ALBERTA AlbertHelp :: localAlbertaFaceNumber[vx];
      // neighborCount is the face number in the actual element
      const int * selfmap  = ALBERTA AlbertHelp :: localAlbertaFaceNumber[nb];

      // copy the three face vertices
      for(int i=0; i<dim; i++)
      {
        const ALBERTA REAL_D & coord = elInfo->coord[selfmap[i]];
        // here consider possible twist
        ALBERTA REAL_D & newcoord    = neighInfo->coord[ neighmap[ facemap[i] ] ];
        for(int j=0; j<dimworld; j++) newcoord[j] = coord[j];
      }
  #endif

      //****************************************
      if (facemap[1] == (facemap[0]+1)%3) {
        return facemap[0];
      }
      // twist
      return facemap[1]-3;
    }
  };
  #endif

  // setup neighbor element with the information of elInfo_
  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>::setupVirtEn() const
  {
    // if this assertion fails then outside was called without checking
    // neighbor first
    assert(neighbor());

    assert( neighborCount_ < dimension+1 );
    // set the neighbor element as element
    // use ALBERTA macro to get neighbour
    neighElInfo_.el = NEIGH(elInfo_->el,elInfo_)[neighborCount_];

    const int vx = elInfo_->opp_vertex[neighborCount_];

    // reset neighbor information
    for( int i = 0; i <= dimension; ++i )
    {
      neighElInfo_.neigh[ i ] = 0;
      neighElInfo_.opp_vertex[ i ] = 127;
    }

    // set origin
    neighElInfo_.neigh[vx] = elInfo_->el;
    neighElInfo_.opp_vertex[vx] = neighborCount_;


  #ifdef CALC_COORD
    // copy the one opposite vertex
    // the same for 2d and 3d
    {
      const ALBERTA REAL_D & coord = elInfo_->opp_coord[neighborCount_];
      ALBERTA REAL_D & newcoord    = neighElInfo_.coord[vx];
      for( int j = 0; j < dimensionworld; ++j )
        newcoord[ j ] = coord[ j ];
    }
  #endif

    // setup coordinates of neighbour elInfo
    twist_ = SetupVirtualNeighbour<GridImp,dimensionworld,dimension>::
             setupNeighInfo(this->grid_,elInfo_,vx,neighborCount_,&neighElInfo_);

    builtNeigh_ = true;
  }

}

#endif
