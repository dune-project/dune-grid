// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_TREEITERATOR_CC
#define DUNE_ALBERTA_TREEITERATOR_CC

#include <dune/grid/albertagrid/treeiterator.cc>

namespace Dune
{

  // AlbertaMarkerVector
  // -------------------

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

  namespace AlbertaMarkerVectorHelp
  {

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


  namespace AlbertaTreeIteratorHelp
  {

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



  // AlbertaGridTreeIterator
  // -----------------------

  //***********************************************************
  //  some template specialization of goNextEntity
  //***********************************************************
  // default implementation, go next elInfo
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALBERTA EL_INFO * AlbertaGridTreeIterator<codim,pitype,GridImp>::
  goNextEntity(ALBERTA TRAVERSE_STACK *stack,ALBERTA EL_INFO *elinfo_old)
  {
    return AlbertaTreeIteratorHelp::GoNextEntity< This, GridImp::dimension, codim >
           ::goNext(*this,stack,elinfo_old);
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
  template< int codim, PartitionIteratorType pitype, class GridImp >
  inline AlbertaGridTreeIterator< codim, pitype, GridImp >
  ::AlbertaGridTreeIterator( const GridImp &grid,
                             int travLevel,
                             int proc,
                             bool leafIt )
    : Base( grid, travLevel, leafIt, true ), // true means end iterator
      level_   (travLevel),
      enLevel_ (travLevel),
      virtualEntity_( this->entityImp() ),
      subEntity_( -1 ),
      vertexMarker_(0),
      okReturn_(false),
      proc_(proc)
  {}


  // Make LevelIterator with point to element from previous iterations
  template< int codim, PartitionIteratorType pitype, class GridImp >
  inline AlbertaGridTreeIterator< codim, pitype, GridImp >
  ::AlbertaGridTreeIterator( const This &other )
    : Base( other.grid_, other.level_, other.leafIt(), (other.vertexMarker_) ? false : true ),
      level_( other.level_ ),
      enLevel_( other.enLevel_ ),
      virtualEntity_( this->entityImp() ),
      manageStack_ (),
      //manageStack_ ( other.manageStack_ ),
      subEntity_( other.subEntity_ ),
      vertexMarker_(other.vertexMarker_),
      okReturn_ (other.okReturn_ ),
      proc_( other.proc_ )
  {
    if(vertexMarker_)
    {
      // if vertexMarker is not NULL then we have a real iterator
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , other.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo, subEntity_ );

      assert( this->grid_.hierarchicIndexSet().index ( *(this->entity_) )
              == this->grid_.hierarchicIndexSet().index ( *(other.entity_) ) );
    }
  }

  // Make LevelIterator with point to element from previous iterations
  template< int codim, PartitionIteratorType pitype, class GridImp >
  inline typename AlbertaGridTreeIterator< codim, pitype, GridImp >::This &
  AlbertaGridTreeIterator<codim,pitype,GridImp>::operator= ( const This &other )
  {
    level_ = other.level_;
    enLevel_ = other.enLevel_;
    //manageStack_ = org.manageStack_;
    subEntity_ =  other.subEntity_;
    vertexMarker_ = (other.vertexMarker_);
    okReturn_ = (other.okReturn_ );

    assert( proc_ == other.proc_ );
    if(vertexMarker_)
    {
      // if vertexMarker is not NULL then we have a real iterator
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , other.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo, subEntity_ );

      assert( this->grid_.hierarchicIndexSet().index ( *(this->entity_) )
              == this->grid_.hierarchicIndexSet().index ( *(other.entity_) ) );
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
      if(!elInfo) return 0;  // if no more Faces, return 0 which leads to end
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
      if(!elInfo) return 0;  // if no more Edges, return
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
      if(!elInfo) return 0;  // if no more Vertices, return
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

    if (stack->stack_used == 0)   /* first call */
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

}

#endif
