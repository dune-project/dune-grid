// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_CC
#define DUNE_ALBERTA_INTERSECTION_CC

#include <dune/grid/albertagrid/intersection.hh>

namespace Dune
{

  template< class GridImp >
  inline AlbertaGridIntersectionIterator<GridImp>::
  AlbertaGridIntersectionIterator(const GridImp & grid,int level) :
    grid_(grid),
    level_ (level),
    neighborCount_ (dim+1),
    elInfo_ (0),
    fakeNeighObj_(LocalGeometryImp()),
    fakeSelfObj_ (LocalGeometryImp()),
    neighGlobObj_(LocalGeometryImp()),
    fakeNeigh_ (grid_.getRealImplementation(fakeNeighObj_) ),
    fakeSelf_  (grid_.getRealImplementation(fakeSelfObj_)  ),
    neighGlob_ (grid_.getRealImplementation(neighGlobObj_)),
    neighElInfo_ () ,
    done_(true)
  {}

  template< class GridImp >
  template< class EntityType >
  inline void AlbertaGridIntersectionIterator<GridImp>::
  first(const EntityType & en ,int level)
  {
    level_ = level;
    neighborCount_ = 0;
    builtNeigh_    = false;
    elInfo_        = en.getElInfo();
    done_          = false;
    this->leafIt_  = en.leafIt();
    assert( elInfo_ );
  }

  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>::
  done ()
  {
    level_ = -1;
    neighborCount_ = dim+1;
    builtNeigh_    = false;
    elInfo_        = 0;
    done_          = true;
  }

  // copy constructor
  template< class GridImp >
  inline AlbertaGridIntersectionIterator<GridImp>::AlbertaGridIntersectionIterator
    (const AlbertaGridIntersectionIterator<GridImp> & org)
    : grid_(org.grid_)
      , level_(org.level_)
      , neighborCount_(org.neighborCount_)
      , builtNeigh_ (false)
      , leafIt_( org.leafIt_ )
      , elInfo_ ( org.elInfo_ )
      , fakeNeighObj_(LocalGeometryImp())
      , fakeSelfObj_ (LocalGeometryImp())
      , neighGlobObj_(LocalGeometryImp())
      , fakeNeigh_ (grid_.getRealImplementation(fakeNeighObj_))
      , fakeSelf_  (grid_.getRealImplementation(fakeSelfObj_ ))
      , neighGlob_ (grid_.getRealImplementation(neighGlobObj_))
      , neighElInfo_()
      , done_ ( org.done_ )
  {}

  // assignment operator
  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>::
  assign (const AlbertaGridIntersectionIterator<GridImp> & org)
  {
    // only assign iterators from the same grid
    assert( &this->grid_ == &(org.grid_));
    level_ =  org.level_;
    neighborCount_ = org.neighborCount_;
    elInfo_ = org.elInfo_;
    builtNeigh_ = false;
    leafIt_ = org.leafIt_;
    done_ = org.done_;
  }


  template< class GridImp >
  inline bool AlbertaGridIntersectionIterator<GridImp>::
  equals (const AlbertaGridIntersectionIterator<GridImp> & i) const
  {
    const ALBERTA EL * e1 = (elInfo_)   ? (elInfo_->el)   : 0;
    const ALBERTA EL * e2 = (i.elInfo_) ? (i.elInfo_->el) : 0;

    return
      ( (e1 == e2)   // check equality of entities which can be both zero
        && (done_ == i.done_)); /// and then also check done status
  }

  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>::increment()
  {
    builtNeigh_ = false;
    // is like go to the next neighbour
    ++neighborCount_;

    // (dim+1) is neigbourCount for end iterators
    if(neighborCount_ > dim)
    {
      this->done();
      return ;
    }

    /*
       // check whether neighbor has same level
       if( (this->neighbor()) && (!this->leafIt()) )
       {
       // if levels are not equal go next
       if( !neighborHasSameLevel () ) increment ();
       }
     */

    return ;
  }

  template< class GridImp >
  inline typename AlbertaGridIntersectionIterator<GridImp>::EntityPointer
  AlbertaGridIntersectionIterator<GridImp>::outside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImpl;

    if(!builtNeigh_)
    {
      assert( elInfo_ );
      // just copy elInfo and then set some values
      std::memcpy(&neighElInfo_,elInfo_,sizeof(ALBERTA EL_INFO));

      setupVirtEn();

      assert( level_ == elInfo_->level );
      assert( (this->leafIt() ) ? (1) : (elInfo_->level == neighElInfo_.level) );
    }

    assert( builtNeigh_ );
    assert( neighElInfo_.el );
    const int neighborLevel = this->grid_.getLevelOfElement( neighElInfo_.el );
    return EntityPointerImpl( this->grid_, neighborLevel, &neighElInfo_, 0 );
  }

  template< class GridImp >
  inline typename AlbertaGridIntersectionIterator<GridImp>::EntityPointer
  AlbertaGridIntersectionIterator<GridImp>::inside () const
  {
    typedef AlbertaGridEntityPointer< 0, GridImp > EntityPointerImpl;
    assert( elInfo_ );
    return EntityPointerImpl( this->grid_, (int)elInfo_->level, elInfo_, 0 );
  }

  template< class GridImp >
  inline int
  AlbertaGridIntersectionIterator<GridImp>::boundaryId () const
  {
    // id of interior intersections is 0
    if( ! boundary() ) return 0;
    assert(elInfo_);
    assert( WALL_BOUND(elInfo_, dim, neighborCount_ ) != 0 );
    return WALL_BOUND(elInfo_, dim, neighborCount_ );
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
  inline const typename AlbertaGridIntersectionIterator<GridImp>::NormalVecType &
  AlbertaGridIntersectionIterator<GridImp>::unitOuterNormal (const LocalCoordType & local) const
  {
    // calculates the outer_normal
    unitNormal_  = this->outerNormal(local);
    unitNormal_ *= (1.0/unitNormal_.two_norm());

    return unitNormal_;
  }

  template<class GridImp>
  inline const typename AlbertaGridIntersectionIterator<GridImp>::NormalVecType &
  AlbertaGridIntersectionIterator<GridImp>::integrationOuterNormal (const LocalCoordType & local) const
  {
    return this->outerNormal(local);
  }


  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator<GridImp>::NormalVecType &
  AlbertaGridIntersectionIterator<GridImp>::outerNormal(const LocalCoordType & local) const
  {
    calcOuterNormal();
    return outNormal_;
  }

  template< class GridImp >
  inline void AlbertaGridIntersectionIterator<GridImp>:: calcOuterNormal() const
  {
    std::cout << "outer_normal() not correctly implemented yet! \n";
    assert(false);
    for(int i=0; i<dimworld; i++)
      outNormal_[i] = 0.0;

    return ;
  }

  template <>
  inline void
  AlbertaGridIntersectionIterator<const AlbertaGrid<2,2> >::calcOuterNormal () const
  {
    assert( elInfo_ );
    const ALBERTA REAL_D & coordOne = grid_.getCoord(elInfo_,(neighborCount_+1)%3);
    const ALBERTA REAL_D & coordTwo = grid_.getCoord(elInfo_,(neighborCount_+2)%3);

    outNormal_[0] = -(coordOne[1] - coordTwo[1]);
    outNormal_[1] =   coordOne[0] - coordTwo[0];
    return ;
  }

  template <>
  inline void AlbertaGridIntersectionIterator<const AlbertaGrid<3,3> >::
  calcOuterNormal () const
  {
    assert( elInfo_ );
    enum { dim = 3 };

    // in this case the orientation is negative, therefore multiply with -1
  #if DIM == 3
    const albertCtype val = (elInfo_->orientation > 0) ? 1.0 : -1.0;
  #else
    const albertCtype val = 1.0;
  #endif

    // neighborCount_ is the local face number
    const int * localFaces = ALBERTA AlbertHelp::localAlbertaFaceNumber[neighborCount_];

    const ALBERTA REAL_D & coordZero = grid_.getCoord(elInfo_, localFaces[0] );
    const ALBERTA REAL_D & coordOne  = grid_.getCoord(elInfo_, localFaces[1] );
    const ALBERTA REAL_D & coordTwo  = grid_.getCoord(elInfo_, localFaces[2] );

    for(int i=0; i<dim; ++i)
    {
      tmpV_[i] = coordOne[i] - coordZero[i];
      tmpU_[i] = coordTwo[i] - coordOne[i];
    }

    // outNormal_ has length 3
    for(int i=0; i<dim; ++i)
      outNormal_[i] = tmpU_[(i+1)%dim] * tmpV_[(i+2)%dim]
                      - tmpU_[(i+2)%dim] * tmpV_[(i+1)%dim];

    outNormal_ *= val;
    return ;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator<GridImp>::LocalGeometry &
  AlbertaGridIntersectionIterator<GridImp>::
  intersectionSelfLocal () const
  {
    fakeSelf_.builtLocalGeom(inside()->geometry(),intersectionGlobal(),
                             elInfo_,neighborCount_);
    return fakeSelfObj_;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator<GridImp>::LocalGeometry &
  AlbertaGridIntersectionIterator<GridImp>::intersectionNeighborLocal () const
  {
    assert(neighbor());

    if(fakeNeigh_.builtLocalGeom(outside()->geometry(),intersectionGlobal(),
                                 &neighElInfo_,neighborCount_)
       )
      return fakeNeighObj_;
    else
    {
      DUNE_THROW(AlbertaError, "intersection_neighbor_local: error occured!");
    }
    return fakeNeighObj_;
  }

  template< class GridImp >
  inline const typename AlbertaGridIntersectionIterator<GridImp>::Geometry &
  AlbertaGridIntersectionIterator<GridImp>::
  intersectionGlobal () const
  {
    assert( elInfo_ );

    if( !neighGlob_.builtGeom( grid_, elInfo_, neighborCount_ ) )
      DUNE_THROW( AlbertaError, "intersectionGlobal: Could not build geometry" );
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

    assert( neighborCount_ < dim+1 );
    // set the neighbor element as element
    // use ALBERTA macro to get neighbour
    neighElInfo_.el = NEIGH(elInfo_->el,elInfo_)[neighborCount_];

    const int vx = elInfo_->opp_vertex[neighborCount_];

    // reset neighbor information
    for(int i=0; i<dim+1; ++i)
    {
      neighElInfo_.neigh[i] = 0;
      neighElInfo_.opp_vertex[i] = 127;
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
      for(int j=0; j<dimworld; j++) newcoord[j] = coord[j];
    }
  #endif

    // setup coordinates of neighbour elInfo
    twist_ = SetupVirtualNeighbour<GridImp,dimworld,dim>::
             setupNeighInfo(this->grid_,elInfo_,vx,neighborCount_,&neighElInfo_);

    builtNeigh_ = true;
  }

}

#endif
