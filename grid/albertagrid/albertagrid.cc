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

namespace Dune
{

  //***************************************************************
  //
  //  --AlbertaGridHierarchicIterator
  //  --HierarchicIterator
  //
  //***************************************************************
  template< class GridImp >
  inline void AlbertaGridHierarchicIterator<GridImp>::
  makeIterator()
  {
    virtualEntity_.setTraverseStack(0);
    virtualEntity_.setElInfo(0,0,0,0);
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const GridImp & grid,
                                int actLevel,
                                int maxLevel)
    : AlbertaGridEntityPointer<0,GridImp> (grid,actLevel,true,true)
      , startLevel_(actLevel)
      , level_ (actLevel)
      , maxlevel_ (maxLevel)
      , virtualEntity_( this->entityImp() )
      , end_ (true)
  {
    makeIterator();
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const GridImp & grid,
                                ALBERTA TRAVERSE_STACK *travStack,int actLevel, int maxLevel, bool leafIt )
    : AlbertaGridEntityPointer<0,GridImp> (grid,actLevel,leafIt,false)
      , startLevel_(actLevel)
      , level_ (actLevel)
      , maxlevel_ ( maxLevel)
      , virtualEntity_( this->entityImp() )
      , end_ (false)
  {
    if(travStack)
    {
      // get new ALBERTA TRAVERSE STACK
      manageStack_.create();

      ALBERTA TRAVERSE_STACK *stack = manageStack_.getStack();

      // cut old traverse stack, kepp only actual element
      ALBERTA copyTraverseStack(stack, travStack );

      // set new traverse level
      if(maxlevel_ < 0)
      {
        std::cout << "WARNING: maxlevel < 0 in AlbertaGridHierarchicIterator! \n";
        // this means, we go until leaf level
        stack->traverse_fill_flag = CALL_LEAF_EL | stack->traverse_fill_flag;
        // exact here has to stand Grid->maxlevel, but is ok anyway
        maxlevel_ = this->grid_.maxLevel();
      }

      // set new traverse level
      stack->traverse_level = maxlevel_;

      virtualEntity_.setTraverseStack(stack);

      // Hier kann ein beliebiges Geometry uebergeben werden,
      // da jedes AlbertGeometry einen Zeiger auf das Macroelement
      // enthaelt.
      ALBERTA EL_INFO * elInfo = firstChild(stack);

      virtualEntity_.setElInfo(elInfo);
    }
    else
    {
      std::cout << "Warning: travStack == NULL in HierarchicIterator(travStack,travLevel) \n";
      makeIterator();
    }
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const AlbertaGridHierarchicIterator<GridImp> & org)
    : AlbertaGridEntityPointer<0,GridImp> (org.grid_,org.level(),true, org.end_ )
      , startLevel_( org.startLevel_ )
      , level_ ( org.level_ )
      , maxlevel_ ( org.maxlevel_ )
      , virtualEntity_( this->entityImp() )
  {
    if( org.virtualEntity_.getElInfo() )
    {
      // get new ALBERTA TRAVERSE STACK
      manageStack_.create();
      ALBERTA TRAVERSE_STACK *stack = manageStack_.getStack();
      // cut old traverse stack, kepp only actual element
      ALBERTA copyTraverseStack(stack, org.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo );
    }
    else
      this->done();
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp> &
  AlbertaGridHierarchicIterator<GridImp>::
  operator = (const AlbertaGridHierarchicIterator<GridImp> & org)
  {
    const_cast<int &> (startLevel_) = org.startLevel_;
    level_    = org.level_;
    maxlevel_ = org.maxlevel_;

    if(org.manageStack_.stackExists())
    {
      // full copy of stack
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      const ALBERTA TRAVERSE_STACK * orgStack = org.manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , orgStack );
    }

    if( org.virtualEntity_.getElInfo() )
      virtualEntity_.setEntity( org.virtualEntity_ );
    else
      this->done();
    return *this;
  }

  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >::increment()
  {
    ALBERTA EL_INFO * nextinfo = recursiveTraverse(manageStack_.getStack());

    if(!nextinfo)
    {
      this->done();
      return;
    }

    virtualEntity_.setElInfo( nextinfo );
    return ;
  }

  template< class GridImp >
  inline ALBERTA EL_INFO *
  AlbertaGridHierarchicIterator<GridImp>::
  firstChild (ALBERTA TRAVERSE_STACK * stack)
  {
    assert(stack);
    assert(stack->elinfo_stack);

    // stack_used ist the actual element
    stack->stack_used = startLevel_+1;

    // info_stack is 0, we want to visit both children
    stack->info_stack[stack->stack_used] = 0;

    ALBERTA EL * el = stack->elinfo_stack[stack->stack_used].el;

    // go down next child
    if(el->child[0] && (stack->traverse_level >
                        (stack->elinfo_stack+stack->stack_used)->level) )
    {
      if(stack->stack_used >= stack->stack_size - 1)
        ALBERTA enlargeTraverseStack(stack);

      int i = stack->info_stack[stack->stack_used];
      el = el->child[i];
      stack->info_stack[stack->stack_used]++;

      // new: go down maxlevel, but fake the elements
      level_++;
      this->grid_.fillElInfo(i, level_, stack->elinfo_stack+stack->stack_used,
                             stack->elinfo_stack+stack->stack_used+1 ,true);

      stack->stack_used++;
      stack->info_stack[stack->stack_used] = 0;
      return (stack->elinfo_stack + stack->stack_used);
    }
    else
    {
      return 0;
    }
  }

  template< class GridImp >
  inline ALBERTA EL_INFO *
  AlbertaGridHierarchicIterator<GridImp>::
  recursiveTraverse(ALBERTA TRAVERSE_STACK * stack)
  {
    // see function
    // static EL_INFO *traverse_leaf_el(TRAVERSE_STACK *stack)
    // Common/traverse_nr_common.cc, line 392
    ALBERTA EL * el=0;

    if(!stack->elinfo_stack)
    {
      /* somethin' wrong */
      return 0;
    }
    else
    {
      // go up until we can go down again
      el = stack->elinfo_stack[stack->stack_used].el;

      // stack->stack_used is actual element in stack
      // stack->info_stack[stack->stack_used] >= 2
      //    means the two children has been visited
      while((stack->stack_used-startLevel_ > 0) &&
            ((stack->info_stack[stack->stack_used] >= 2)
             || (el->child[0] == 0)
             || ( stack->traverse_level <=
                  (stack->elinfo_stack+stack->stack_used)->level)) )
      {
        stack->stack_used--;
        el = stack->elinfo_stack[stack->stack_used].el;
        level_ = stack->elinfo_stack[stack->stack_used].level;
      }

      // goto next father is done by other iterator and not our problem
      if(stack->stack_used-startLevel_ < 1)
      {
        return 0;
      }
    }

    // go down next child
    if(el->child[0] && (stack->traverse_level >
                        (stack->elinfo_stack+stack->stack_used)->level) )
    {
      if(stack->stack_used >= stack->stack_size - 1)
        ALBERTA enlargeTraverseStack(stack);

      int i = stack->info_stack[stack->stack_used];
      el = el->child[i];
      stack->info_stack[stack->stack_used]++;

      // new: go down maxlevel, but fake the elements
      level_++;
      this->grid_.fillElInfo(i, level_, stack->elinfo_stack+stack->stack_used,
                             stack->elinfo_stack+stack->stack_used+1 ,true);

      stack->stack_used++;
      stack->info_stack[stack->stack_used] = 0;
    }
    else
    {
      return 0;
    }

    return (stack->elinfo_stack + stack->stack_used);
  } // recursive traverse over all childs
    // end AlbertaGridHierarchicIterator

  //***************************************************************
  //
  //  --AlbertaGridIntersectionIterator
  //  --IntersectionIterator
  //
  //***************************************************************

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
      ( (e1 == e2) // check equality of entities which can be both zero
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
  // end IntersectionIterator


  //*******************************************************
  //
  // --AlbertaGridTreeIterator
  // --TreeIterator
  // --LevelIterator
  //
  //*******************************************************

  namespace AlbertaTreeIteratorHelp {

    // for elements
    template <class IteratorImp, int dim>
    struct GoNextEntity<IteratorImp,dim,0>
    {
      static inline ALBERTA EL_INFO *
      goNext(IteratorImp & it, ALBERTA TRAVERSE_STACK *stack , ALBERTA EL_INFO *elinfo_old)
      {
        return it.goNextElInfo(stack,elinfo_old);
      }
    };

    // for faces
    template <class IteratorImp, int dim>
    struct GoNextEntity<IteratorImp,dim,1>
    {
      static inline ALBERTA EL_INFO *
      goNext(IteratorImp & it, ALBERTA TRAVERSE_STACK *stack , ALBERTA EL_INFO *elinfo_old)
      {
        return it.goNextFace(stack,elinfo_old);
      }
    };

    // for vertices
    template <class IteratorImp, int dim>
    struct GoNextEntity<IteratorImp,dim,dim>
    {
      static inline ALBERTA EL_INFO *
      goNext(IteratorImp & it, ALBERTA TRAVERSE_STACK *stack , ALBERTA EL_INFO *elinfo_old)
      {
        return it.goNextVertex(stack,elinfo_old);
      }
    };

    // for edges in 3d
    template <class IteratorImp>
    struct GoNextEntity<IteratorImp,3,2>
    {
      static inline ALBERTA EL_INFO *
      goNext(IteratorImp & it, ALBERTA TRAVERSE_STACK *stack , ALBERTA EL_INFO *elinfo_old)
      {
        return it.goNextEdge(stack,elinfo_old);
      }
    };
  } // end namespace AlbertaTreeIteratorHelp


  //***********************************************************
  //  some template specialization of goNextEntity
  //***********************************************************
  // default implementation, go next elInfo
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextEntity(ALBERTA TRAVERSE_STACK *stack,ALBERTA EL_INFO *elinfo_old)
  {
    return AlbertaTreeIteratorHelp :: GoNextEntity<ThisType,GridImp::dimension,codim>::
           goNext(*this,stack,elinfo_old);
  }
  //***************************************

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline void AlbertaGridTreeIterator<codim,pitype,GridImp>::
  makeIterator()
  {
    level_ = 0;
    enLevel_ = 0;
    subEntity_ = -1;
    vertexMarker_ = 0;

    virtualEntity_.setTraverseStack(0);
    virtualEntity_.setElInfo( 0, 0 );
  }

  // Make LevelIterator with point to element from previous iterations
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline AlbertaGridTreeIterator<codim,pitype,GridImp>::
  AlbertaGridTreeIterator(const GridImp & grid, int travLevel,int proc, bool leafIt )
    : AlbertaGridEntityPointer<codim,GridImp> (grid,travLevel,leafIt,true) // true means end iterator
      , level_   (travLevel)
      , enLevel_ (travLevel)
      , virtualEntity_( this->entityImp() )
      , subEntity_( -1 )
      , vertexMarker_(0)
      , okReturn_(false)
      , proc_(proc)
  {}

  // Make LevelIterator with point to element from previous iterations
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline AlbertaGridTreeIterator<codim,pitype,GridImp>::
  AlbertaGridTreeIterator(const AlbertaGridTreeIterator<codim,pitype,GridImp> & org)
    : AlbertaGridEntityPointer<codim,GridImp> (org.grid_,org.level_, org.leafIt() , (org.vertexMarker_) ? false : true)
      , level_   (org.level_)
      , enLevel_ (org.enLevel_)
      , virtualEntity_( this->entityImp() )
      , manageStack_ ()
      //, manageStack_ ( org.manageStack_ )
      , subEntity_( org.subEntity_ )
      , vertexMarker_(org.vertexMarker_)
      , okReturn_ (org.okReturn_ )
      , proc_(org.proc_)
  {
    if(vertexMarker_)
    {
      // if vertexMarker is not NULL then we have a real iterator
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , org.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo, subEntity_ );

      assert( this->grid_.hierarchicIndexSet().index ( *(this->entity_) )
              == this->grid_.hierarchicIndexSet().index ( *(org.entity_) ) );
    }
  }

  // Make LevelIterator with point to element from previous iterations
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline AlbertaGridTreeIterator<codim,pitype,GridImp> &
  AlbertaGridTreeIterator<codim,pitype,GridImp>::operator =
    (const AlbertaGridTreeIterator<codim,pitype,GridImp> & org)
  {
    level_ = org.level_;
    enLevel_ = org.enLevel_;
    //manageStack_ = org.manageStack_;
    subEntity_ =  org.subEntity_;
    vertexMarker_ = (org.vertexMarker_);
    okReturn_ = (org.okReturn_ );

    assert( proc_ == org.proc_ );
    if(vertexMarker_)
    {
      // if vertexMarker is not NULL then we have a real iterator
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , org.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo, subEntity_ );

      assert( this->grid_.hierarchicIndexSet().index ( *(this->entity_) )
              == this->grid_.hierarchicIndexSet().index ( *(org.entity_) ) );
    }

    return *this;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline AlbertaGridTreeIterator<codim,pitype,GridImp>::
  AlbertaGridTreeIterator(const GridImp & grid,
                          const AlbertaMarkerVector * vertexMark,
                          int travLevel, int proc, bool leafIt)
    : AlbertaGridEntityPointer<codim,GridImp> (grid,travLevel,leafIt,false)
      , level_ (travLevel) , enLevel_(travLevel)
      , virtualEntity_( this->entityImp() )
      , subEntity_( -1 )
      , vertexMarker_(0)
      , okReturn_ (false)
      , proc_(proc)
  {
    ALBERTA MESH * mesh = this->grid_.getMesh();

    if( mesh && ((travLevel >= 0) && (travLevel <= this->grid_.maxLevel())) )
    {
      vertexMarker_ = vertexMark;

      ALBERTA FLAGS travFlags = FILL_ALL(mesh); //FILL_COORDS | FILL_NEIGH;

      // CALL_LEAF_EL is not used anymore
      travFlags = travFlags | CALL_LEAF_EL_LEVEL;

      // get traverse_stack
      manageStack_.create();

      virtualEntity_.setTraverseStack(manageStack_.getStack());

      // diese Methode muss neu geschrieben werden, da man
      // die ParentElement explizit speichern moechte.
      ALBERTA EL_INFO* elInfo =
        goFirstElement(manageStack_.getStack(), mesh, travLevel,travFlags);

      virtualEntity_.setElInfo(elInfo, subEntity_ );
    }
    else
    {
      // create empty iterator
      makeIterator();
    }
  }

  // gehe zum naechsten Element, wie auch immer
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline void AlbertaGridTreeIterator<codim,pitype,GridImp>::increment()
  {
    ALBERTA EL_INFO * nextinfo = goNextEntity(manageStack_.getStack(),virtualEntity_.getElInfo());

    if(!nextinfo)
    {
      this->done();
      return ;
    }

    virtualEntity_.setElInfo( nextinfo , subEntity_ );

    return ;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextFace(ALBERTA TRAVERSE_STACK *stack, ALBERTA EL_INFO *elInfo)
  {
    // go next Element, if face_ > numberOfVertices, then go to next elInfo
    // face_ is set to -1 by constructor
    ++subEntity_;
    if( subEntity_ >= (dim+1)) // dim+1 Faces
    {
      // we have checked all faces on this element,
      // therefore go to next element
      elInfo = goNextElInfo(stack, elInfo);
      subEntity_ = 0;
      if(!elInfo) return 0; // if no more Faces, return 0 which leads to end
    }

    // check elInfo pointer before we start anything
    assert(elInfo);

    if( !this->leafIt() )
    {
      // go next, if Vertex is not treated on this Element
      if(vertexMarker_->faceNotOnElement(
           this->grid_.getElementNumber(elInfo->el),
           this->grid_.getFaceNumber(elInfo->el,subEntity_)))
      {
        elInfo = goNextFace(stack,elInfo);
      }
    }
    else
    {
      // get neighbour of this element
      const ALBERTA EL * neighbour = NEIGH(elInfo->el,elInfo)[ subEntity_ ];
      if( neighbour )
      {
        // get element
        const ALBERTA EL * el = elInfo->el;
        assert( el );

        // if true we must go to next face
        // when element number is small then go next because now the face is
        // reached on the element with the largest number
        bool goWeida = this->grid_.getElementNumber( el ) < this->grid_.getElementNumber( neighbour ) ;

        // additional check for level iterators
        if(goWeida && ! this->leafIt() )
        {
          // if we should go weida then only go
          // if neighbours are on the same level
          goWeida = (this->grid_.getLevelOfElement( neighbour ) == level_);
        }

        // the face was reached before therefore go to next face
        if(goWeida) elInfo = goNextFace(stack,elInfo);
      }
    }

    return elInfo;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextEdge(ALBERTA TRAVERSE_STACK *stack, ALBERTA EL_INFO *elInfo)
  {
    // go next Element, Edge 0
    // treat Edge like Faces
    // edge_ is set to -1 by constructor
    ++subEntity_;
    if( subEntity_ >= 6 ) // in 3d only 6 Edges
    {
      elInfo = goNextElInfo(stack, elInfo);
      subEntity_ = 0;
      if(!elInfo) return 0; // if no more Edges, return
    }

    assert(elInfo);

    // go next, if Vertex is not treated on this Element
    if(vertexMarker_->edgeNotOnElement(
         this->grid_.getElementNumber(elInfo->el),
         this->grid_.getEdgeNumber(elInfo->el,subEntity_)))
    {
      elInfo = goNextEdge(stack,elInfo);
    }

    return elInfo;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextVertex(ALBERTA TRAVERSE_STACK *stack, ALBERTA EL_INFO *elInfo)
  {
    // go next Element, Vertex 0
    // treat Vertices like Faces
    // vertex_ is set to -1 by constructor
    ++subEntity_;
    if( subEntity_ >= (dim+1)) // dim+1 Vertices
    {
      elInfo = goNextElInfo(stack, elInfo);
      subEntity_ = 0;
      if(!elInfo) return 0; // if no more Vertices, return
    }

    assert(elInfo);

    // go next, if Vertex is not treated on this Element
    if(vertexMarker_->vertexNotOnElement(
         this->grid_.getElementNumber(elInfo->el),
         this->grid_.getVertexNumber(elInfo->el,subEntity_)))
    {
      elInfo = goNextVertex(stack,elInfo);
    }

    return elInfo;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goFirstElement(ALBERTA TRAVERSE_STACK *stack,ALBERTA MESH *mesh, int level,
                 ALBERTA FLAGS fill_flag)
  {
    FUNCNAME("goFirstElement");

    if (!stack)
    {
      ALBERTA_ERROR("no traverse stack\n");
      return 0;
    }

    stack->traverse_mesh      = mesh;
    stack->traverse_level     = level;
    stack->traverse_fill_flag = fill_flag;

    if (stack->stack_size < 1)
      enlargeTraverseStack(stack);

    for (int i=0; i<stack->stack_size; i++)
      stack->elinfo_stack[i].fill_flag = fill_flag & FILL_ALL(mesh);

    stack->elinfo_stack[0].mesh = stack->elinfo_stack[1].mesh = mesh;

    if (fill_flag & CALL_LEAF_EL_LEVEL)
    {
      ALBERTA_TEST_EXIT(level >= 0) ("invalid level: %d\n",level);
    }

    stack->traverse_mel = 0;
    stack->stack_used   = 0;
    stack->el_count     = 0;

    // go to first enInfo, therefore goNextElInfo
    ALBERTA EL_INFO * elinfo = goNextElInfo(stack,0);
    // for codim 0 we are done at this point
    if( codim == 0 ) return elinfo;

    // if iterate over faces, edges or vertices
    // we have to go the normal way, because
    // is not preset that the first face lies on the first element
    // therefore face_,edge_,vertex_ == -1 at the start
    return goNextEntity(stack,elinfo);
  }


  // --travNext
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextElInfo(ALBERTA TRAVERSE_STACK *stack, ALBERTA EL_INFO *elinfo_old)
  {
    assert( stack );
    //assert( elinfo_old );
    assert( (stack->stack_used) ?
            (elinfo_old == stack->elinfo_stack+stack->stack_used) :
            (elinfo_old == 0) );
    {
      // overloaded traverse_leaf_el_level, is not implemened in ALBERTA yet
      ALBERTA EL_INFO  * elinfo = traverseElLevel(stack);

      // if leafIt_ == false go to elements only on desired level
      if((elinfo) && (! this->leafIt()) )
      {
        if(elinfo->level == stack->traverse_level)
          okReturn_ = true;

        while(!okReturn_)
        {
          elinfo = traverseElLevel(stack);
          if(!elinfo) okReturn_ = true;
        }
        ++(stack->el_count);
      }

      // set new level for Entity
      if((elinfo) && (this->leafIt()) ) enLevel_ = elinfo->level;

      return(elinfo);
    }
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  traverseElLevel(ALBERTA TRAVERSE_STACK *stack)
  {
    FUNCNAME("traverseElLevel");
    ALBERTA EL *el;
    int i;
    okReturn_ = false;

    if (stack->stack_used == 0) /* first call */
    {
      stack->traverse_mel = FIRST_MACRO_EL(stack->traverse_mesh);
      if (stack->traverse_mel == nil) return(nil);

      stack->stack_used = 1;

      ALBERTA fillMacroInfo(stack, stack->traverse_mel,
                            stack->elinfo_stack+stack->stack_used,level_);

      stack->info_stack[stack->stack_used] = 0;

      el = stack->elinfo_stack[stack->stack_used].el;
      if ((el == nil) || (el->child[0] == nil)) {
        return(stack->elinfo_stack+stack->stack_used);
      }
    }
    else
    {
      el = stack->elinfo_stack[stack->stack_used].el;

      /* go up in tree until we can go down again */
      while((stack->stack_used > 0) &&
            ((stack->info_stack[stack->stack_used] >= 2) || (el->child[0]==nil)
             || ( stack->traverse_level <=
                  (stack->elinfo_stack+stack->stack_used)->level)) )
      {
        stack->stack_used--;
        el = stack->elinfo_stack[stack->stack_used].el;
      }
      /* goto next macro element */
      if (stack->stack_used < 1)
      {
        stack->traverse_mel = NEXT_MACRO_EL(stack->traverse_mesh,stack->traverse_mel);
        if (stack->traverse_mel == nil) return(nil);

        stack->stack_used = 1;

        ALBERTA fillMacroInfo(stack, stack->traverse_mel,
                              stack->elinfo_stack+stack->stack_used,level_);

        stack->info_stack[stack->stack_used] = 0;

        el = stack->elinfo_stack[stack->stack_used].el;
        if ((el == nil) || (el->child[0] == nil))
        {
          return(stack->elinfo_stack+stack->stack_used);
        }
      }
    }

    /* go down tree until leaf oder level*/
    while (el->child[0] &&
           ( stack->traverse_level > (stack->elinfo_stack+stack->stack_used)->level))
    {
      if(stack->stack_used >= stack->stack_size-1)
        enlargeTraverseStack(stack);

      i = stack->info_stack[stack->stack_used];
      el = el->child[i];
      stack->info_stack[stack->stack_used]++;

      this->grid_.fillElInfo(i, level_, stack->elinfo_stack+stack->stack_used,
                             stack->elinfo_stack+stack->stack_used+1, false , this->leafIt() );

      stack->stack_used++;

      ALBERTA_TEST_EXIT(stack->stack_used < stack->stack_size)
        ("stack_size=%d too small, level=(%d,%d)\n",
        stack->stack_size, stack->elinfo_stack[stack->stack_used].level);

      stack->info_stack[stack->stack_used] = 0;

      if(stack->traverse_level == (stack->elinfo_stack+stack->stack_used)->level)
      {
        okReturn_ = true;
      }
    }

    return(stack->elinfo_stack+stack->stack_used);
  }

  //*************************************************************************
  //  end AlbertaGridTreeIterator
  //*************************************************************************

  template <int dim, class GridImp>
  inline AlbertaGridHierarchicIterator<GridImp>
  AlbertaGridEntity <0,dim,GridImp>::hbegin(int maxlevel) const
  {
    // Kopiere alle Eintraege des stack, da man im Stack weiterlaeuft und
    // sich deshalb die Werte anedern koennen, der elinfo_stack bleibt jedoch
    // der gleiche, deshalb kann man auch nur nach unten, d.h. zu den Kindern
    // laufen
    return AlbertaGridHierarchicIterator<GridImp> (grid_,travStack_,level(),maxlevel, this->leafIt() );
  }

  template <int dim, class GridImp>
  inline AlbertaGridHierarchicIterator<GridImp>
  AlbertaGridEntity <0,dim,GridImp>::hend(int maxlevel) const
  {
    AlbertaGridHierarchicIterator<GridImp> it(grid_,level(),maxlevel);
    return it;
  }

  template <int dim, class GridImp>
  inline typename AlbertaGridEntity <0,dim,GridImp>::AlbertaGridLeafIntersectionIteratorType
  AlbertaGridEntity <0,dim,GridImp>::ileafbegin() const
  {
#ifndef NDEBUG
    for(int i=0; i<GridImp::dimension+1; ++i)
    {
      if(elInfo_->opp_vertex[i] == 127 )
      {
        DUNE_THROW(NotImplemented,"AlbertaGridIntersectionIterator::first: do not create IntersectionIterators on outside entities, not implemented yet!");
      }
    }
#endif
    return AlbertaGridLeafIntersectionIteratorType(grid_,*this,level(),false);
  }

  template <int dim, class GridImp>
  inline typename AlbertaGridEntity <0,dim,GridImp>::AlbertaGridLeafIntersectionIteratorType
  AlbertaGridEntity <0,dim,GridImp>::ileafend() const
  {
    return AlbertaGridLeafIntersectionIteratorType(grid_,*this,level(),true);
  }

  //*********************************************************************
  //
  //  AlbertMarkerVertex
  //
  //*********************************************************************
  inline bool AlbertaMarkerVector::
  vertexNotOnElement(const int elIndex, const int vertex) const
  {
    return (vec_[ vertex ] != elIndex);
  }

  inline bool AlbertaMarkerVector::
  edgeNotOnElement(const int elIndex, const int edge) const
  {
    return (edgevec_[ edge ] != elIndex);
  }

  inline bool AlbertaMarkerVector::
  faceNotOnElement(const int elIndex, const int face) const
  {
    assert( facevec_.size () > 0 );
    return (facevec_[ face ] != elIndex);
  }

  namespace AlbertaMarkerVectorHelp {

    // only for 3d calc edges markers
    template <class GridType>
    struct MarkFaces
    {
      template <class ArrayType>
      inline static void mark(GridType & grid , ArrayType & vec, const ALBERTA EL * el, int elindex)
      {
        enum { dim = GridType :: dimension };
        // we have dim+1 faces
        for(int i=0; i<dim+1; ++i)
        {
          // the 1 is the difference between dim=3 and codim=2
          int num = grid.getFaceNumber(el ,i);
          if( vec[num] == -1 ) vec[num] = elindex;
        }
      }
    };

    template <class GridType, int dim>
    struct MarkEdges
    {
      template <class ArrayType>
      inline static void mark(GridType & grid , ArrayType & vec, const ALBERTA EL * el, int elindex)
      {}
    };

    // only for 3d calc edges markers
    template <class GridType>
    struct MarkEdges<GridType,3>
    {
      template <class ArrayType>
      inline static void mark(GridType & grid , ArrayType & vec, const ALBERTA EL * el, int elindex)
      {
        // in 3d 6 edges
        for(int i=0; i<6; ++i)
        {
          // the 1 is the difference between dim=3 and codim=2
          int num = grid.getEdgeNumber(el ,i);
          if( vec[num] == -1 ) vec[num] = elindex;
        }
      }
    };

  } // end namespace AlbertaMarkerVectorHelp

  template <class GridType>
  inline void AlbertaMarkerVector::markNewVertices(GridType &grid, int level)
  {
    assert( meLevel_ == true );
    enum { dim      = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: HierarchicIndexSet HIndexSet;

    const HIndexSet & hset = grid.hierarchicIndexSet();
    int nvx = hset.size(dim);
    int fce = hset.size(1);

    {
      ArrayType & vec     = vec_;
      if((int) vec.size() < nvx) vec.resize( nvx + vxBufferSize_ );
      const int vecSize = vec.size();
      for(int i=0; i<vecSize; i++) vec[i] = -1;

      ArrayType & facevec     = facevec_;
      if((int) facevec.size() < fce) facevec.resize( fce + vxBufferSize_ );
      const int facevecSize = facevec.size();
      for(int i=0; i<facevecSize; i++) facevec[i] = -1;

      ArrayType & edgevec = edgevec_;
      if(dim > 2)
      {
        int edg = hset.size(dim-1);
        if((int) edgevec.size() < edg) edgevec.resize( edg + vxBufferSize_ );
        const int edgevecSize = edgevec.size();
        for(int i=0; i<edgevecSize; i++) edgevec[i] = -1;
      }

      typedef typename GridType::template Codim<0>::LevelIterator LevelIteratorType;
      LevelIteratorType endit = grid.template lend<0> (level);
      for(LevelIteratorType it = grid.template lbegin<0> (level); it != endit; ++it)
      {
        const ALBERTA EL * el =
          (grid.getRealImplementation(*it)).getElInfo()->el;

        int elindex = grid.getElementNumber(el);
        for(int local=0; local<dim+1; local++)
        {
          int num = grid.getVertexNumber(el,local); // vertex num
          if( vec[num] == -1 ) vec[num] = elindex;
        }

        // mark edges for this element
        // in 3d 6 edges
        AlbertaMarkerVectorHelp::MarkFaces<GridType>::mark(grid,facevec, el,elindex );

        // mark edges for this element
        // in 3d 6 edges
        AlbertaMarkerVectorHelp::MarkEdges<GridType,dim>::mark(grid,edgevec, el,elindex );
      }
      // remember the number of entity on level and codim = 0
    }
    up2Date_ = true;
  }

  // mark vertices and edges using leaf iterator
  template <class GridType>
  inline void AlbertaMarkerVector::markNewLeafVertices(GridType &grid)
  {
    assert( meLevel_ == false );
    enum { dim      = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    int nvx = grid.hierarchicIndexSet().size(dim);

    {
      ArrayType & vec = vec_;
      if((int) vec.size() < nvx) vec.resize( nvx + vxBufferSize_ );

      // the edge marking is only needed in 3d
      ArrayType & edgevec = edgevec_;
      if( dim > 2 )
      {
        int edg = grid.hierarchicIndexSet().size(dim-1);
        if((int) edgevec.size() < edg) edgevec.resize( edg + vxBufferSize_ );
        const int edgevecSize = edgevec.size();
        for(int i=0; i<edgevecSize; i++) edgevec[i] = -1;
      }

      const int vecSize = vec.size();
      for(int i=0; i<vecSize; i++) vec[i] = -1;

      typedef typename GridType::template Codim<0>::LeafIterator IteratorType;
      IteratorType endit = grid.template leafend<0> ();
      for(IteratorType it = grid.template leafbegin<0> (); it != endit; ++it)
      {
        const ALBERTA EL * el =
          (grid.getRealImplementation(*it)).getElInfo()->el;

        int elindex = grid.hierarchicIndexSet().index(*it);
        for(int local=0; local<dim+1; local++)
        {
          int num = el->dof[local][0]; // vertex num
          if( vec[num] == -1 ) vec[num] = elindex;
        }

        // mark edges for this element
        AlbertaMarkerVectorHelp::MarkEdges<GridType,dim>::mark(grid,edgevec,el,elindex);
      }
      // remember the number of entity on leaf level and codim = 0
    }
    up2Date_ = true;
  }

  inline void AlbertaMarkerVector::print() const
  {
    {
      if(vec_.size() > 0)
      {
        std::cout << "\nEntries " << vec_.size() << std::endl;
        const int size = vec_.size();
        for(int i=0; i<size; ++i)
        {
          std::cout << "Vx " << i << " visited on Element " << vec_[i] << std::endl;
        }
      }
    }
  }

  //***********************************************************************
  //
  // --AlbertaGrid
  // --Grid
  //
  //***********************************************************************
  template < int dim, int dimworld >
  inline AlbertaGrid < dim, dimworld >::AlbertaGrid()
    : mesh_ (0)
      , comm_()
      , maxlevel_ (0) , wasChanged_ (false)
      , vertexMarkerLeaf_(false) // creates LeafMarkerVector
      , nv_ (dim+1) , dof_ (0)
      , myRank_ (0)
      , hIndexSet_(*this)
      , globalIdSet_(*this)
      , levelIndexVec_(MAXL,0)
      , leafIndexSet_ (0)
      , geomTypes_()
      , sizeCache_ (0)
      , coarsenMarked_(0)
      , refineMarked_(0)
      , lockPostAdapt_(false)
  {
    // create vector with geom types
    makeGeomTypes();

    for(int i=0; i<AlbertHelp::numOfElNumVec; i++) dofvecs_.elNumbers[i] = 0;
    dofvecs_.elNewCheck = 0;
#ifndef CALC_COORD
    dofvecs_.coords     = 0;
#endif
  }

  template< int dim, int dimworld >
  inline void AlbertaGrid< dim, dimworld > :: makeGeomTypes ()
  {
    // we have dim+1 codims
    geomTypes_.resize( dim+1 );

    // for each codim create geom type
    for( int codim = 0; codim <= dim; ++codim )
    {
      const GeometryType type( GeometryType :: simplex, dim - codim );
      geomTypes_[ codim ].push_back( type );
    }
  }

  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::initGrid(int proc)
  {
    ALBERTA AlbertHelp::getDofVecs<dimworld> (&dofvecs_);

    calcExtras();

    wasChanged_ = true;

    macroVertices_.resize( mesh_->n_vertices );

    LeafDataType::initLeafDataValues(mesh_,proc);

    calcExtras();
  }

  template < int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld >
  :: AlbertaGrid( const std::string macroTriangFilename )
    : mesh_( 0 ),
      comm_(),
      maxlevel_( 0 ),
      wasChanged_( false ),
      vertexMarkerLeaf_( false ), // creates LeafMarkerVector
      nv_( dim+1 ),
      dof_( 0 ),
      myRank_( -1 ),
      hIndexSet_( *this ),
      globalIdSet_( *this ),
      levelIndexVec_( MAXL, 0 ),
      leafIndexSet_ ( 0 ),
      geomTypes_(),
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

    const char * MacroTriangFilename = macroTriangFilename.c_str();
    assert( MacroTriangFilename );

    bool makeNew = true;
    {
      std::fstream file (MacroTriangFilename,std::ios::in);
      if(!file) DUNE_THROW(AlbertaIOError,"could not open grid file " << MacroTriangFilename);

      std::basic_string <char> str,str1;
      file >> str1; str = str1.assign(str1,0,3);
      // With that Albert MacroTriang starts DIM or DIM_OF_WORLD
      if (str != "DIM") makeNew = false;
      file.close();
    }

    ALBERTA AlbertHelp::initIndexManager_elmem_cc(indexStack_);

    if(makeNew)
    {
      ALBERTA AlbertHelp :: initBndStack( &bndStack_ );
      // create mesh
      mesh_ = ALBERTA AlbertHelp::createMesh<dim> ("AlbertaGrid", MacroTriangFilename);

      ALBERTA AlbertHelp :: removeBndStack ();

      initGrid(0);
    }
    else
    {
      derr<<"Couldn't read grid file '"<<macroTriangFilename<<"' because it's not in ALBERTA macro triangulation format! \n";
      DUNE_THROW(NotImplemented,"Constructor reading backup file not implemented!");
    }
    std::cout << "AlbertaGrid<"<<dim<<","<<dimworld<<"> created from macro grid file '" << macroTriangFilename << "'. \n\n";
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

    for(int i=0; i<AlbertHelp::numOfElNumVec; i++)
    {
      if(dofvecs_.elNumbers[i]) ALBERTA free_dof_int_vec(dofvecs_.elNumbers[i]);
      dofvecs_.elNumbers[i] = 0;
    }

    if(dofvecs_.elNewCheck) ALBERTA free_dof_int_vec(dofvecs_.elNewCheck);dofvecs_.elNewCheck = 0;
#ifndef CALC_COORD
    if(dofvecs_.coords ) ALBERTA free_dof_real_d_vec(dofvecs_.coords);dofvecs_.coords = 0;
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
    if(mesh_) ALBERTA free_mesh(mesh_);mesh_ = 0;

    // delete all created boundary structures
    while ( !bndStack_.empty() )
    {
      ALBERTA BOUNDARY * obj = bndStack_.top();
      bndStack_.pop();
      if( obj ) delete obj;
    }
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
  template <class GridImp, class EntityProvider, int dim , int codim >
  struct GetNewEntity
  {
    typedef typename SelectEntityImp<codim,dim,GridImp>::EntityObject EntityObject;
    typedef typename SelectEntityImp<codim,dim,GridImp>::EntityImp EntityImp;
    static EntityObject *
    getNewEntity(GridImp & grid, EntityProvider &enp , int level, bool leafIt )
    {
      return new EntityObject (EntityImp(grid,level,leafIt));
    }

    static void freeEntity (EntityProvider &enp , EntityObject * en)
    {
      if(en) delete en;
    }
  };

  // specialisation for codim 0 uses stack
  template <class GridImp, class EntityProvider, int dim>
  struct GetNewEntity<GridImp,EntityProvider,dim,0>
  {
    typedef typename SelectEntityImp<0,dim,GridImp>::EntityObject EntityObject;
    typedef typename SelectEntityImp<0,dim,GridImp>::EntityImp EntityImp;
    static EntityObject *
    getNewEntity(GridImp & grid, EntityProvider &enp , int level, bool leafIt )
    {
      // return object from stack
      return enp.getNewObjectEntity(grid,(EntityImp *)0,level,leafIt);
    }

    static void freeEntity (EntityProvider &enp , EntityObject * en)
    {
      enp.freeObjectEntity(en);
    }
  };

  template < int dim   , int dimworld >
  template < int codim >
  inline typename SelectEntityImp<codim,dim,const AlbertaGrid<dim,dimworld> >::EntityObject *
  AlbertaGrid < dim, dimworld >::getNewEntity (int level, bool leafIt ) const
  {
    return GetNewEntity<const MyType, EntityProvider, dim, codim > :: getNewEntity(*this,entityProvider_,level,leafIt);
  }

  template < int dim   , int dimworld >
  template < int codim >
  inline void AlbertaGrid < dim, dimworld >::
  freeEntity (typename SelectEntityImp<codim,dim,const MyType>::EntityObject * en) const
  {
    GetNewEntity<const MyType, EntityProvider, dim, codim > :: freeEntity(entityProvider_,en);
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

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::preAdapt()
  {
    return (coarsenMarked_ > 0);
  }

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::postAdapt()
  {
    assert( (leafIndexSet_) ? (mesh_->n_elements == leafIndexSet_->size(0) ?   1 : 0) : 1);
    assert( (leafIndexSet_) ? (mesh_->n_vertices == leafIndexSet_->size(dim) ? 1 : 0) : 1);
#if DIM == 3
    //assert( (leafIndexSet_ && dim == 3) ? (mesh_->n_edges == leafIndexSet_->size(dim-1) ?  1 :0) :1);
    assert( (leafIndexSet_ && dim == 3) ? (mesh_->n_faces == leafIndexSet_->size(1) ? 1 : 0) : 1);
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
    ALBERTA AlbertHelp::set2positive(dofvecs_.elNewCheck);

    return wasChanged_;
  }

  template<int dim, int dimworld>
  inline bool AlbertaGrid < dim, dimworld >::
  mark( int refCount , const typename Traits::template Codim<0>::EntityPointer & ep ) const
  {
    return this->mark(refCount,*ep);
  }

  //--mark
  template<int dim, int dimworld>
  inline bool AlbertaGrid < dim, dimworld >::
  mark( int refCount , const typename Traits::template Codim<0>::Entity & ep ) const
  {
    // if not leaf entity, leaf method
    if( !ep.isLeaf() ) return false;

    ALBERTA EL_INFO * elInfo = (this->getRealImplementation(ep)).getElInfo();
    if(!elInfo) return false;
    ALBERTA EL * element = elInfo->el;
    assert( element );

    // mark for refinement
    if( refCount > 0)
    {
      element->mark = refCount;
      int factor = 2;
      for(int i=0; i<refCount; ++i) factor *= 2;
      refineMarked_ += factor;
      return true;
    }

    // mark for coarsening
    if( refCount < 0)
    {
      element->mark = refCount;
      ++coarsenMarked_;
      return true;
    }

    // mark for none
    element->mark = 0;
    return true;
  }

  // --getMark
  template<int dim, int dimworld>
  inline int AlbertaGrid < dim, dimworld >::
  getMark( const typename Traits::template Codim<0>::EntityPointer & ep ) const
  {
    return this->getMark( *ep );
  }

  template<int dim, int dimworld>
  inline int AlbertaGrid < dim, dimworld >::
  getMark( const typename Traits::template Codim<0>::Entity & en ) const
  {
    const ALBERTA EL_INFO * elInfo = (this->getRealImplementation(en)).getElInfo();
    assert( elInfo );
    assert( elInfo->el );
    return elInfo->el->mark;
  }

  // --adapt
  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::adapt()
  {
    unsigned char flag;
    bool refined = false;
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
    ALBERTA AlbertHelp::set2positive ( dofvecs_.elNewCheck );

    flag = ALBERTA AlbertRefine ( mesh_ );
    refined = (flag == 0) ? false : true;

    if(preAdapt()) // true if a least on element is marked for coarseing
      flag = ALBERTA AlbertCoarsen( mesh_ );

    if(!refined)
    {
      wasChanged_ = (flag == 0) ? false : true;
    }
    else
      wasChanged_ = true;

    if(wasChanged_)
    {
      calcExtras();
    }

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

    ALBERTA AlbertHelp::AdaptRestrictProlongHandler < MyType , COType >
    handler ( *this,
              father, this->getRealImplementation(father) ,
              son   , this->getRealImplementation(son) ,
              tmprpop );

    ALBERTA AlbertHelp::MeshCallBack & callBack = ALBERTA
                                                  AlbertHelp::MeshCallBack::instance();

    callBack.setPointers(mesh_,handler);

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

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::checkElNew (const ALBERTA EL *el) const
  {
    // if element is new then entry in dofVec is 1
    return (elNewVec_[el->dof[dof_][nv_]] < 0);
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::maxLevel() const
  {
    return maxlevel_;
  }

  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::global_size (int codim) const
  {
    if(codim == dim) return mesh_->n_vertices;
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
    hIndexSet_.updatePointers(dofvecs_);

#ifndef CALC_COORD
    coordsVec_ = (dofvecs_.coords)->vec;      assert(coordsVec_);
#endif
    elNewVec_  = (dofvecs_.elNewCheck)->vec;  assert(elNewVec_);
    elAdmin_   = dofvecs_.elNumbers[0]->fe_space->admin;

    // see Albert Doc. , should stay the same
    const_cast<int &> (nv_)  = elAdmin_->n0_dof[CENTER];
    const_cast<int &> (dof_) = elAdmin_->mesh->node[CENTER];
  }


  template < int dim, int dimworld >
  inline int AlbertaGrid < dim, dimworld >::getLevelOfElement (const ALBERTA EL *el) const
  {
    assert( el );
    // return the elements level which is the absolute value of the entry
    return std::abs( elNewVec_ [el->dof[dof_][nv_]] );
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
    maxlevel_ = ALBERTA AlbertHelp::calcMaxAbsoluteValueOfVector( dofvecs_.elNewCheck );
    assert( maxlevel_ >= 0);
    assert( maxlevel_ < MAXL);

#ifndef NDEBUG
    int mlvl = ALBERTA AlbertHelp::calcMaxLevel(mesh_,dofvecs_.elNewCheck);
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

  template < int dim, int dimworld >  template <GrapeIOFileFormatType ftype>
  inline bool AlbertaGrid < dim, dimworld >::
  writeGrid (const std::basic_string<char> filename, albertCtype time ) const
  {
    switch(ftype)
    {
    case xdr   : return writeGridXdr (filename , time );
    case ascii :
    {
      // write leaf grid as macro triangulation
      int ret = ALBERTA write_macro( mesh_ , filename.c_str() );
      return (ret == 1) ? true : false;
    }
    default :
    {
      DUNE_THROW(AlbertaError,"wrong FileType in AlbertaGrid::writeGrid!");
      return false;
    }
    }
  }

  template < int dim, int dimworld >  template <GrapeIOFileFormatType ftype>
  inline bool AlbertaGrid < dim, dimworld >::
  readGrid (const std::basic_string<char> filename, albertCtype &time )
  {
    switch(ftype)
    {
    case xdr   : return readGridXdr   (filename , time );
    case ascii : return readGridAscii (filename , time );
    default : {
      DUNE_THROW(AlbertaError,"wrong FileType in AlbertaGrid::readGrid!");
      return false;
    }
    }
  }

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::
  writeGridXdr (const std::basic_string<char> filename, albertCtype time ) const
  {
    std::string ownerfile(filename);
    if(filename.size() > 0)
    {
      ownerfile += "_own";
    }
    else
      DUNE_THROW(AlbertaIOError, "no filename given in writeGridXdr ");

    // strore element numbering to file
    for(int i=0; i<AlbertHelp::numOfElNumVec; i++)
    {
      std::string elnumfile(filename);
      elnumfile += "_num_c";
      char tmpchar[16]; sprintf(tmpchar,"%d",i);
      elnumfile += tmpchar;
      ALBERTA write_dof_int_vec_xdr(dofvecs_.elNumbers[i],elnumfile.c_str());
    }

    const char * fn = filename.c_str();
    int flag = ALBERTA write_mesh_xdr (mesh_ , fn , time);
    return (flag == 1) ? true : false;
  }

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::
  readGridXdr (const std::basic_string<char> filename, albertCtype & time )
  {
    // remove all old stuff
    // to be reivised
    //removeMesh();

    const char * fn = filename.c_str();

    ALBERTA AlbertHelp :: initBndStack( &bndStack_ );
    mesh_ = (ALBERTA read_mesh_xdr (fn , &time ,
                                    LeafDataType::initLeafData,
                                    ALBERTA AlbertHelp::initBoundary) );
    ALBERTA AlbertHelp :: removeBndStack ();

    if (mesh_ == 0)
      DUNE_THROW(AlbertaIOError, "could not open grid file " << filename);

    // read element numbering from file
    std::string ownerfile (filename);
    if(filename.size() > 0)
    {
      ownerfile += "_own";
    }
    else
      return false;

    for(int i=0; i<AlbertHelp::numOfElNumVec; i++)
    {
      std::string elnumfile( filename );
      char tmpchar[16]; sprintf(tmpchar,"%d",i);
      elnumfile += "_num_c"; elnumfile += tmpchar;
      const char * elnumfn = elnumfile.c_str();
      dofvecs_.elNumbers[i] = ALBERTA read_dof_int_vec_xdr(elnumfn, mesh_ , 0 );
    }

    // make the rest of the dofvecs
    ALBERTA AlbertHelp::makeTheRest<dimworld>(&dofvecs_);

    // restore level information for each element by traversing the mesh
    ALBERTA AlbertHelp::restoreElNewCheck( mesh_ , dofvecs_.elNewCheck );

    // make vectors know in grid and hSet
    arrangeDofVec();

    // calc maxlevel and indexOnLevel and so on
    calcExtras();

    // set el_index of index manager to max element index
    for(int i=0; i<ALBERTA AlbertHelp::numOfElNumVec; i++)
    {
      int maxIdx = ALBERTA AlbertHelp::calcMaxIndex( dofvecs_.elNumbers[i] );
      indexStack_[i].setMaxIndex(maxIdx);
    }

    return true;
  }

  template < int dim, int dimworld >
  inline bool AlbertaGrid < dim, dimworld >::readGridAscii
    (const std::basic_string<char> filename, albertCtype & time )
  {
    removeMesh(); // delete all objects

    ALBERTA AlbertHelp :: initBndStack( &bndStack_ );
    // create mesh
    mesh_ = ALBERTA AlbertHelp::createMesh<dim> ("AlbertaGrid", filename.c_str());

    ALBERTA AlbertHelp :: removeBndStack ();

    time = 0.0;

    // unset up2Dat status, if lbegin is called then this status is updated
    for(int l=0; l<MAXL; l++) vertexMarkerLevel_[l].unsetUp2Date();

    // unset up2Dat status, if leafbegin is called then this status is updated
    vertexMarkerLeaf_.unsetUp2Date();

    ALBERTA AlbertHelp::initIndexManager_elmem_cc(indexStack_);

    initGrid(myRank_);
    return true;
  }

  // if defined some debugging test were made that reduce the performance
  // so they were switch off normaly

  //#define DEBUG_FILLELINFO
  //*********************************************************************
  //  fillElInfo 2D
  //*********************************************************************
  typedef U_CHAR ALBERTA_CHAR;

  template<int dim, int dimworld>
  inline void AlbertaGrid<dim,dimworld >::
  firstNeigh(const int ichild, const ALBERTA EL_INFO *elinfo_old,
             ALBERTA EL_INFO *elinfo, const bool leafLevel) const
  {
#ifdef CALC_COORD
    // old stuff
    const ALBERTA REAL_D * old_opp_coord  = elinfo_old->opp_coord;
    const ALBERTA REAL_D * old_coord      = elinfo_old->coord;

    // new stuff
    ALBERTA REAL_D * opp_coord = elinfo->opp_coord;
#endif
    ALBERTA ALBERTA_CHAR * opp_vertex     = elinfo->opp_vertex;
    ALBERTA EL ** neigh = NEIGH(el,elinfo);

    assert(neigh != 0);
    assert(neigh == NEIGH(el,elinfo));

    const int onechi = 1-ichild;

    // first nb 0 of new elinfo
    {
      ALBERTA EL * nb = NEIGH(elinfo_old->el,elinfo_old)[2];
      if(nb)
      {
        // if NULL then nonconforme refinement
        assert(nb->child[0] != 0);
        ALBERTA EL * nextNb = nb->child[onechi];
        opp_vertex[ichild] = onechi;
#ifdef CALC_COORD
        for(int i=0; i<dimworld; i++)
          opp_coord[ichild][i]  = old_opp_coord[2][i];
#endif
        // if we have children, we could go down
        if(nextNb->child[0])
        {
          // if level of neighbour to small, do down once more
          // but only one level down because of conformity
          if( leafLevel )
          {
            nextNb = nextNb->child[ichild];
            opp_vertex[ichild] = 2;
#ifdef CALC_COORD
            for(int i=0; i<dimworld; i++)
              opp_coord[ichild][i] += old_coord[ichild][i];
            for(int i=0; i<dimworld; i++)
              opp_coord[ichild][i] *= 0.5;
#endif
          }
        }
        neigh[ichild] = nextNb;
      }
      else
      {
        // if no neighbour then children have no neighbour
        neigh[ichild] = 0;
      }
    }
  }

  template<int dim, int dimworld>
  inline void AlbertaGrid<dim,dimworld >::
  secondNeigh(const int ichild, const ALBERTA EL_INFO *elinfo_old,
              ALBERTA EL_INFO *elinfo, const bool leafLevel) const
  {
    // old stuff
#ifdef CALC_COORD
    const ALBERTA REAL_D * old_coord      = elinfo_old->coord;

    // new stuff
    ALBERTA REAL_D * opp_coord = elinfo->opp_coord;
#endif
    ALBERTA ALBERTA_CHAR * opp_vertex     = elinfo->opp_vertex;
    ALBERTA EL ** neigh = NEIGH(el,elinfo);

    assert(neigh != 0);
    assert(neigh == NEIGH(el,elinfo));

    const int onechi = 1-ichild;
    // nb 1 of new elinfo, always the same
    {
      ALBERTA EL * nextNb = elinfo_old->el->child[onechi];
      opp_vertex[onechi] = ichild;
#ifdef CALC_COORD
      for(int i=0; i<dimworld; i++)
        opp_coord[onechi][i]  = old_coord[onechi][i];
#endif
      // check if children exists
      if(nextNb->child[0] )
      {
        if( leafLevel )
        {
          nextNb = nextNb->child[onechi];
          opp_vertex[onechi] = 2;
#ifdef CALC_COORD
          for(int i=0; i<dimworld; i++)
            opp_coord[onechi][i] += old_coord[2][i];
          for(int i=0; i<dimworld; i++)
            opp_coord[onechi][i] *= 0.5;
#endif
        }
      }
      neigh[onechi] = nextNb;
    }
  }

  template<int dim, int dimworld>
  inline void AlbertaGrid<dim,dimworld >::
  thirdNeigh(const int ichild, const ALBERTA EL_INFO *elinfo_old,
             ALBERTA EL_INFO *elinfo, const bool leafLevel) const
  {
    // old stuff
    const ALBERTA ALBERTA_CHAR * old_opp_vertex = elinfo_old->opp_vertex;
#ifdef CALC_COORD
    const ALBERTA REAL_D * old_opp_coord  = elinfo_old->opp_coord;
    const ALBERTA REAL_D * old_coord      = elinfo_old->coord;

    // new stuff
    ALBERTA REAL_D * opp_coord = elinfo->opp_coord;
#endif
    ALBERTA ALBERTA_CHAR * opp_vertex     = elinfo->opp_vertex;
    ALBERTA EL ** neigh = NEIGH(el,elinfo);

    assert(neigh != 0);
    assert(neigh == NEIGH(el,elinfo));

    const int onechi = 1-ichild;
    // nb 2 of new elinfo
    {
      ALBERTA EL * nb = NEIGH(elinfo_old->el,elinfo_old)[onechi];
      if(nb)
      {
        const int vx = old_opp_vertex[onechi];
        opp_vertex[2] = vx;
#ifdef CALC_COORD
        for(int i=0; i<dimworld; i++)
          opp_coord[2][i] = old_opp_coord[onechi][i];
#endif
        if((vx == 2) || (nb->child[0] == 0))
        {
          // means the neighbour has the same refinement edge like our child
          neigh[2] = nb;
        }
        else
        {
          assert(nb->child[0] != 0);
          neigh[2] = nb->child[1-vx];
          opp_vertex[2] = 2;
#ifdef CALC_COORD
          const int vxind = (vx == ichild) ? ichild : 2;
          for(int i=0; i<dimworld; i++)
            opp_coord[2][i] += old_coord[vxind][i];
          for(int i=0; i<dimworld; i++)
            opp_coord[2][i] *= 0.5;
#endif
        }
      }
      else
      {
        // if no neighbour then children have no neighbour
        neigh[2] = 0;
      }
    }

  }

  template<int dim, int dimworld>
  inline void AlbertaGrid<dim,dimworld >::
  fillElInfo(int ichild, int actLevel , const ALBERTA EL_INFO *elinfo_old,
             ALBERTA EL_INFO *elinfo, bool hierarchical, bool leaf) const
  {
#if DUNE_ALBERTA_VERSION >= 0x200
    // the alberta version of filling an EL_INFO structure
    ALBERTA fill_elinfo(ichild,elinfo_old,elinfo);
#else

    // old stuff
    ALBERTA EL * el_old = elinfo_old->el;
    assert(el_old != 0);
    assert(el_old->child[0] != 0);

    // new stuff
    // set new element
    ALBERTA EL * el  = el_old->child[ichild];
    elinfo->el = el;

    ALBERTA FLAGS fill_flag = elinfo_old->fill_flag;

    elinfo->macro_el  = elinfo_old->macro_el;
    elinfo->fill_flag = fill_flag;
    elinfo->mesh      = elinfo_old->mesh;
    elinfo->parent    = el_old;
    elinfo->level     = elinfo_old->level + 1;

    // calculate the coordinates
    if (fill_flag & FILL_COORDS)
    {
      if (el_old->new_coord)
      {
        for (int j = 0; j < dimworld; j++)
          elinfo->coord[2][j] = el_old->new_coord[j];
      }
      else
      {
        for (int j = 0; j < dimworld; j++)
          elinfo->coord[2][j] =
            0.5 * (elinfo_old->coord[0][j] + elinfo_old->coord[1][j]);
      }

      // set the other coord
      for (int j = 0; j < dimworld; j++)
      {
        elinfo->coord[ichild  ][j] = elinfo_old->coord[2][j];
        elinfo->coord[1-ichild][j] = elinfo_old->coord[ichild][j];
      }
    }

    // calculate neighbours
    // neighbor info does not work for hierarchical walk
    if(fill_flag & FILL_NEIGH && ! hierarchical )
    {
      // allow to go down on neighbour more than once
      // if the following condition is satisfied
      //const bool leafLevel = ((el->child[0] == 0) || (elinfo->level < actLevel));
      const bool leafLevel = (leaf) ? true : ((el->child[0] == 0) && (elinfo->level < actLevel));

      firstNeigh (ichild,elinfo_old,elinfo,leafLevel);
      secondNeigh(ichild,elinfo_old,elinfo,leafLevel);
      thirdNeigh (ichild,elinfo_old,elinfo,leafLevel);
    }

    // boundary calculation
    if (fill_flag & FILL_BOUND)
    {
      if (elinfo_old->boundary[2])
        elinfo->bound[2] = elinfo_old->boundary[2]->bound;
      else
        elinfo->bound[2] = INTERIOR;

      elinfo->bound[ichild]   = elinfo_old->bound[2];
      elinfo->bound[1-ichild] = elinfo_old->bound[ichild];
      elinfo->boundary[ichild] = elinfo_old->boundary[2];
      elinfo->boundary[1-ichild] = nil;
      elinfo->boundary[2] = elinfo_old->boundary[1-ichild];
    }

#endif
  } // end Grid::fillElInfo 2D


  //***********************************************************************
  // fillElInfo 3D
  //***********************************************************************
#if DIM == 3
  template <>
  inline void AlbertaGrid<3,3>::
  fillElInfo(int ichild, int actLevel , const ALBERTA EL_INFO *elinfo_old,
             ALBERTA EL_INFO *elinfo, bool hierarchical, bool leaf) const
  {
    ALBERTA fill_elinfo(ichild,elinfo_old,elinfo);
  }
#endif

  template < int dim, int dimworld >
  inline void AlbertaGrid < dim, dimworld >::setNewCoords
    (const FieldVector<albertCtype, dimworld> & trans, const albertCtype scalar)
  {
    static FieldVector<albertCtype, dimworld> trans_(0.0);
    static albertCtype scalar_ (1.0);

    for(int i=0; i<macroVertices_.size(); i++)
      macroVertices_[i] = 0;

    for(ALBERTA MACRO_EL * mel = FIRST_MACRO_EL(mesh_); mel; mel = NEXT_MACRO_EL(mesh_,mel) )
    {
      for(int i=0; i<dim+1; i++)
      {
        int dof = mel->el->dof[i][0];
        // visit each coord only once
        if( macroVertices_[dof] != 1)
        {
          macroVertices_[dof] = 1;
          for(int j=0; j<dimworld; j++)
          {
            mel->coord[i][j] -= trans_[j];
            mel->coord[i][j] /= scalar_;

            mel->coord[i][j] *= scalar;
            mel->coord[i][j] += trans[j];
          }
        }
      }
    }

    for (int i=0; i<dimworld; i++)
      trans_[i] = trans[i];

    scalar_ = scalar;
  }

} // namespace Dune

#undef ALBERTA_CHAR
#endif
