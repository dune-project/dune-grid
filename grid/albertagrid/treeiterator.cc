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

    template< class Grid, int codim >
    class MarkSubEntities
    {
      typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;

      static const int dimension = Grid::dimension;
      static const int numSubEntities = Alberta::NumSubEntities< dimension, codim >::value;

      typedef Alberta::ElementInfo< dimension > ElementInfo;

    public:
      template< class Array >
      static void apply ( const Grid &grid, Array &vec, const ElementInfo &elementInfo )
      {
        const HierarchicIndexSet &hIndexSet = grid.hierarchicIndexSet();

        const int index = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
        for( int i = 0; i < numSubEntities; ++i )
        {
          const int subIndex = hIndexSet.template subIndex< codim >( elementInfo, i );
          if( vec[ subIndex ] < 0 )
            vec[ subIndex ] = index;
        }
      }
    };

  } // end namespace AlbertaMarkerVectorHelp


  template< class Grid >
  inline void AlbertaMarkerVector::markNewVertices ( const Grid &grid, int level )
  {
    typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;
    typedef typename Grid::template Codim< 0 >::LevelIterator LevelIterator;

    const int dim = Grid::dimension;

    assert( meLevel_ == true );

    const HierarchicIndexSet &hIndexSet = grid.hierarchicIndexSet();

    size_t nvx = hIndexSet.size( dim );
    int fce = hIndexSet.size( 1 );

    std::vector< int > &vec = vec_;
    if( vec.size() < nvx )
      vec.resize( nvx + vxBufferSize_ );

    const int vecSize = vec.size();
    for(int i=0; i<vecSize; i++) vec[i] = -1;

    std::vector< int > &facevec = facevec_;
    if((int) facevec.size() < fce) facevec.resize( fce + vxBufferSize_ );
    const int facevecSize = facevec.size();
    for(int i=0; i<facevecSize; i++) facevec[i] = -1;

    std::vector< int > &edgevec = edgevec_;
    if( dim > 2 )
    {
      int edg = hIndexSet.size(dim-1);
      if((int) edgevec.size() < edg) edgevec.resize( edg + vxBufferSize_ );
      const int edgevecSize = edgevec.size();
      for(int i=0; i<edgevecSize; i++) edgevec[i] = -1;
    }

    const LevelIterator endit = grid.template lend< 0 >( level );
    for( LevelIterator it = grid.template lbegin< 0 >( level ); it != endit; ++it )
    {
      const Alberta::ElementInfo< dim > &elementInfo
        = Grid::getRealImplementation( *it ).elementInfo();
      AlbertaMarkerVectorHelp::MarkSubEntities< Grid, dim >::apply( grid, vec, elementInfo );
      AlbertaMarkerVectorHelp::MarkSubEntities< Grid, 1 >::apply( grid, facevec, elementInfo );
      if( dim == 3 )
        AlbertaMarkerVectorHelp::MarkSubEntities< Grid, 2 >::apply( grid, edgevec, elementInfo );
    }

    up2Date_ = true;
  }


  // mark vertices and edges using leaf iterator
  template< class Grid >
  inline void AlbertaMarkerVector::markNewLeafVertices( const Grid &grid )
  {
    typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;
    typedef typename Grid::template Codim< 0 >::LeafIterator LeafIterator;

    const int dim = Grid::dimension;

    assert( meLevel_ == false );

    int nvx = grid.hierarchicIndexSet().size(dim);

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

    const LeafIterator endit = grid.template leafend< 0 >();
    for( LeafIterator it = grid.template leafbegin< 0 >(); it != endit; ++it )
    {
      const Alberta::ElementInfo< dim > &elementInfo
        = Grid::getRealImplementation( *it ).elementInfo();
      AlbertaMarkerVectorHelp::MarkSubEntities< Grid, dim >::apply( grid, vec, elementInfo );
      if( dim == 3 )
        AlbertaMarkerVectorHelp::MarkSubEntities< Grid, 2 >::apply( grid, edgevec, elementInfo );
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
    : Base( grid ),
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
    // it is ok to set the invalid ElementInfo
    entityImp().setElement( elementInfo, subEntity_ );
  }


  // Make LevelIterator with point to element from previous iterations
  template< int codim, class GridImp, bool leafIterator >
  inline AlbertaGridTreeIterator< codim, GridImp, leafIterator >
  ::AlbertaGridTreeIterator ( const GridImp &grid,
                              int travLevel )
    : Base( grid ),
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
    // it is ok to set the invalid ElementInfo
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
    return (leafIterator ? elementInfo.isLeaf() : (level_ == elementInfo.level()));
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

    const HierarchicIndexSet &hIndexSet = grid().hierarchicIndexSet();
    if( leafIterator )
    {
      const ALBERTA EL *neighbor = elementInfo.elInfo().neigh[ subEntity_ ];
      if( neighbor != NULL )
      {
        // face is reached from element with largest number
        const int elIndex = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
        const int nbIndex = hIndexSet.template subIndex< 0 >( neighbor, 0 );
        if( elIndex < nbIndex )
          goNextFace( elementInfo );
      }
    }
    else
    {
      const int elIndex = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
      const int faceIndex = hIndexSet.template subIndex< 1 >( elementInfo, subEntity_ );
      assert( vertexMarker_ != 0 );
      if( vertexMarker_->faceNotOnElement( elIndex, faceIndex ) )
        goNextFace( elementInfo );
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

    const HierarchicIndexSet &hIndexSet = grid().hierarchicIndexSet();
    const int elIndex = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
    const int edgeIndex = hIndexSet.template subIndex< 2 >( elementInfo, subEntity_ );
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

    const HierarchicIndexSet &hIndexSet = grid().hierarchicIndexSet();
    const int elIndex = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
    const int vertexIndex = hIndexSet.template subIndex< dimension >( elementInfo, subEntity_ );
    assert( vertexMarker_ != 0 );
    if( vertexMarker_->vertexNotOnElement( elIndex, vertexIndex ) )
      goNextVertex( elementInfo );
  }

}

#endif
