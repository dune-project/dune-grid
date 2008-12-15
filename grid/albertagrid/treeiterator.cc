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
    assert( facevec_.size() > 0 );
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
      std::vector< int > &vec = vec_;
      if((int) vec.size() < nvx) vec.resize( nvx + vxBufferSize_ );
      const int vecSize = vec.size();
      for(int i=0; i<vecSize; i++) vec[i] = -1;

      std::vector< int > &facevec = facevec_;
      if((int) facevec.size() < fce) facevec.resize( fce + vxBufferSize_ );
      const int facevecSize = facevec.size();
      for(int i=0; i<facevecSize; i++) facevec[i] = -1;

      std::vector< int > &edgevec = edgevec_;
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
      std::vector< int > &vec = vec_;
      if((int) vec.size() < nvx) vec.resize( nvx + vxBufferSize_ );

      // the edge marking is only needed in 3d
      std::vector< int > &edgevec = edgevec_;
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



  // AlbertaTreeIterratorHelp
  // ------------------------

  namespace AlbertaTreeIteratorHelp
  {

    // for elements
    template< class IteratorImp, int dim >
    struct GoNextEntity< IteratorImp, dim, 0 >
    {
      typedef typename IteratorImp::ElementInfo ElementInfo;

      static void goNext ( IteratorImp &it, ElementInfo &elementInfo )
      {
        it.goNextElement( elementInfo );
      }
    };

    // for faces
    template <class IteratorImp, int dim>
    struct GoNextEntity<IteratorImp,dim,1>
    {
      typedef typename IteratorImp::ElementInfo ElementInfo;

      static void goNext ( IteratorImp &it, ElementInfo &elementInfo )
      {
        it.goNextFace( elementInfo );
      }
    };

    // for vertices
    template <class IteratorImp, int dim>
    struct GoNextEntity<IteratorImp,dim,dim>
    {
      typedef typename IteratorImp::ElementInfo ElementInfo;

      static void goNext ( IteratorImp &it, ElementInfo &elementInfo )
      {
        it.goNextVertex( elementInfo );
      }
    };

    // for edges in 3d
    template <class IteratorImp>
    struct GoNextEntity<IteratorImp,3,2>
    {
      typedef typename IteratorImp::ElementInfo ElementInfo;

      static void goNext ( IteratorImp &it, ElementInfo &elementInfo )
      {
        it.goNextEdge( elementInfo );
      }
    };

  } // end namespace AlbertaTreeIteratorHelp



  // AlbertaGridTreeIterator
  // -----------------------

  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNextEntity ( ElementInfo &elementInfo )
  {
    return AlbertaTreeIteratorHelp::GoNextEntity< This, GridImp::dimension, codim >
           ::goNext( *this, elementInfo );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >::makeIterator ()
  {
    level_ = 0;
    subEntity_ = -1;
    vertexMarker_ = 0;

    entityImp().clearElement();
  }


  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ( const GridImp &grid,
                              const AlbertaMarkerVector *vertexMark,
                              int travLevel )
    : Base( grid, travLevel, leafIterator, false ),
      level_( travLevel ),
      subEntity_( (codim == 0 ? 0 : -1) ),
      macroIterator_( grid.meshPointer().begin() ),
      vertexMarker_( vertexMark )
  {
    ElementInfo elementInfo = *macroIterator_;
    if( codim == 0 )
      nextElementStop( elementInfo );
    else
      goNextEntity( elementInfo );
    if( !elementInfo )
      this->done();
    else
      entityImp().setElement( elementInfo, subEntity_ );
  }


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ( const GridImp &grid,
                              int travLevel )
    : Base( grid, travLevel, leafIterator, true ), // true means end iterator
      level_( travLevel ),
      subEntity_( -1 ),
      macroIterator_( grid.meshPointer().end() ),
      vertexMarker_( 0 )
  {}


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator( const This &other )
    : Base( other ),
      level_( other.level_ ),
      subEntity_( other.subEntity_ ),
      macroIterator_( other.macroIterator_ ),
      vertexMarker_( other.vertexMarker_ )
  {}


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline typename AlbertaGridTreeIterator< codim, GridImp, leafIterator >::This &
  AlbertaGridTreeIterator< codim, GridImp, leafIterator >::operator= ( const This &other )
  {
    Base::operator=( other );

    level_ = other.level_;
    subEntity_ =  other.subEntity_;
    macroIterator_ = other.macroIterator_;
    vertexMarker_ = other.vertexMarker_;

    return *this;
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >::increment ()
  {
    ElementInfo elementInfo = entityImp().elementInfo_;
    goNextEntity ( elementInfo );
    if( !elementInfo )
      this->done();
    else
      entityImp().setElement( elementInfo, subEntity_ );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::nextElement ( ElementInfo &elementInfo )
  {
    if( elementInfo.isLeaf() )
    {
      while( (elementInfo.level() > 0) && (elementInfo.indexInFather() == 1) )
        elementInfo = elementInfo.father();
      if( elementInfo.level() == 0 )
      {
        ++macroIterator_;
        elementInfo = *macroIterator_;
      }
      else
        elementInfo = elementInfo.father().child( 1 );
    }
    else
      elementInfo = elementInfo.child( 0 );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::nextElementStop ( ElementInfo &elementInfo )
  {
    while( !(!elementInfo || stopAtElement( elementInfo )) )
      nextElement( elementInfo );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline bool AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::stopAtElement ( const ElementInfo &elementInfo )
  {
    if( !elementInfo )
      return true;
    return (this->leafIt() ? elementInfo.isLeaf() : (level_ == elementInfo.level()));
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNextElement ( ElementInfo &elementInfo )
  {
    nextElement( elementInfo );
    nextElementStop( elementInfo );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNextFace ( ElementInfo &elementInfo )
  {
    ++subEntity_;
    if( subEntity_ >= numSubEntities )
    {
      subEntity_ = 0;
      nextElement( elementInfo );
      nextElementStop( elementInfo );
      if( !elementInfo )
        return;
    }

    ALBERTA EL *el = elementInfo.el();
    assert( el );
    if( !this->leafIt() )
    {
      const int elIndex = this->grid_.getElementNumber( el );
      const int faceIndex = this->grid_.getFaceNumber( el, subEntity_ );
      assert( vertexMarker_ != 0 );
      if( vertexMarker_->faceNotOnElement( elIndex, faceIndex ) )
        goNextFace( elementInfo );
    }
    else
    {
      // get neighbour of this element
      const ALBERTA EL *neighbor = elementInfo.elInfo().neigh[ subEntity_ ];
      if( neighbor != NULL )
      {
        const int elIndex = this->grid_.getElementNumber( el );
        const int nbIndex = this->grid_.getElementNumber( neighbor );

        // when element number is small then go next because now the face is
        // reached on the element with the largest number
        if( elIndex < nbIndex )
          goNextFace( elementInfo );
      }
    }
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNextEdge ( ElementInfo &elementInfo )
  {
    ++subEntity_;
    if( subEntity_ >= numSubEntities )
    {
      subEntity_ = 0;
      nextElement( elementInfo );
      nextElementStop( elementInfo );
      if( !elementInfo )
        return;
    }

    ALBERTA EL *el = elementInfo.el();
    assert( el );
    const int elIndex = this->grid_.getElementNumber( el );
    const int edgeIndex = this->grid_.getEdgeNumber( el, subEntity_ );
    assert( vertexMarker_ != 0 );
    if( vertexMarker_->edgeNotOnElement( elIndex, edgeIndex ) )
      goNextEdge( elementInfo );
  }


  template< int codim, class GridImp, bool leafIterator >
  inline void AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::goNextVertex ( ElementInfo &elementInfo )
  {
    ++subEntity_;
    if( subEntity_ >= numSubEntities )
    {
      subEntity_ = 0;
      nextElement( elementInfo );
      nextElementStop( elementInfo );
      if( !elementInfo )
        return;
    }

    ALBERTA EL *el = elementInfo.el();
    assert( el );
    const int elIndex = this->grid_.getElementNumber( el );
    const int vertexIndex = this->grid_.getVertexNumber( el, subEntity_ );
    assert( vertexMarker_ != 0 );
    if( vertexMarker_->vertexNotOnElement( elIndex, vertexIndex ) )
      goNextVertex( elementInfo );
  }

}

#endif
