// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_TREEITERATOR_CC
#define DUNE_ALBERTA_TREEITERATOR_CC

#include <dune/grid/albertagrid/treeiterator.cc>

namespace Dune
{

  // AlbertaMarkerVector
  // -------------------

  template< int dim, int dimworld >
  inline bool AlbertaMarkerVector< dim, dimworld >
  ::vertexNotOnElement(const int elIndex, const int vertex) const
  {
    const int codim = dimension;
    assert( marker_[ codim ].size() > 0 );
    return (marker_[ codim ][ vertex ] != elIndex);
  }

  template< int dim, int dimworld >
  inline bool AlbertaMarkerVector< dim, dimworld >
  ::edgeNotOnElement(const int elIndex, const int edge) const
  {
    const int codim = 2;
    assert( marker_[ codim ].size() > 0 );
    return (marker_[ codim ][ edge ] != elIndex);
  }

  template< int dim, int dimworld >
  inline bool AlbertaMarkerVector< dim, dimworld >
  ::faceNotOnElement(const int elIndex, const int face) const
  {
    const int codim = 1;
    assert( marker_[ codim ].size() > 0 );
    return (marker_[ codim ][ face ] != elIndex);
  }



  template< int dim, int dimworld >
  template< int codim >
  class AlbertaMarkerVector< dim, dimworld >::MarkSubEntities
  {
    static const int numSubEntities = Alberta::NumSubEntities< dimension, codim >::value;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  public:
    template< class Array >
    static void apply ( const HierarchicIndexSet &hIndexSet,
                        Array (&marker)[ dimension + 1],
                        const ElementInfo &elementInfo )
    {
      Array &vec = marker[ codim ];

      const int index = hIndexSet.template subIndex< 0 >( elementInfo, 0 );
      for( int i = 0; i < numSubEntities; ++i )
      {
        const int subIndex = hIndexSet.template subIndex< codim >( elementInfo, i );
        if( vec[ subIndex ] < 0 )
          vec[ subIndex ] = index;
      }
    }
  };


  template< int dim, int dimworld >
  inline void AlbertaMarkerVector< dim, dimworld >
  ::markNewVertices ( const Grid &grid, int level )
  {
    typedef typename Grid::template Codim< 0 >::LevelIterator LevelIterator;

    assert( meLevel_ == true );

    const HierarchicIndexSet &hIndexSet = grid.hierarchicIndexSet();

    size_t nvx = hIndexSet.size( dimension );
    int fce = hIndexSet.size( 1 );

    std::vector< int > &vec = marker_[ dimension ];
    if( vec.size() < nvx )
      vec.resize( nvx + vxBufferSize_ );

    const int vecSize = vec.size();
    for(int i=0; i<vecSize; i++) vec[i] = -1;

    std::vector< int > &facevec = marker_[ 1 ];
    if((int) facevec.size() < fce) facevec.resize( fce + vxBufferSize_ );
    const int facevecSize = facevec.size();
    for(int i=0; i<facevecSize; i++) facevec[i] = -1;

    std::vector< int > &edgevec = marker_[ 2 ];
    if( dimension > 2 )
    {
      int edg = hIndexSet.size( dimension - 1 );
      if((int) edgevec.size() < edg) edgevec.resize( edg + vxBufferSize_ );
      const int edgevecSize = edgevec.size();
      for(int i=0; i<edgevecSize; i++) edgevec[i] = -1;
    }

    const LevelIterator endit = grid.template lend< 0 >( level );
    for( LevelIterator it = grid.template lbegin< 0 >( level ); it != endit; ++it )
    {
      const Alberta::ElementInfo< dimension > &elementInfo
        = Grid::getRealImplementation( *it ).elementInfo();
      Alberta::ForLoop< MarkSubEntities, 1, dimension >::apply( hIndexSet, marker_, elementInfo );
    }

    up2Date_ = true;
  }


  // mark vertices and edges using leaf iterator
  template< int dim, int dimworld >
  inline void AlbertaMarkerVector< dim, dimworld >
  ::markNewLeafVertices ( const Grid &grid )
  {
    typedef typename Grid::template Codim< 0 >::LeafIterator LeafIterator;

    assert( meLevel_ == false );

    const HierarchicIndexSet &hIndexSet = grid.hierarchicIndexSet();
    int nvx = hIndexSet.size( dimension );

    std::vector< int > &vec = marker_[ dimension ];
    if((int) vec.size() < nvx) vec.resize( nvx + vxBufferSize_ );

    // the edge marking is only needed in 3d
    std::vector< int > &edgevec = marker_[ 2 ];
    if( dimension > 2 )
    {
      int edg = hIndexSet.size( dimension-1 );
      if((int) edgevec.size() < edg) edgevec.resize( edg + vxBufferSize_ );
      const int edgevecSize = edgevec.size();
      for(int i=0; i<edgevecSize; i++) edgevec[i] = -1;
    }

    const int vecSize = vec.size();
    for(int i=0; i<vecSize; i++) vec[i] = -1;

    const LeafIterator endit = grid.template leafend< 0 >();
    for( LeafIterator it = grid.template leafbegin< 0 >(); it != endit; ++it )
    {
      const Alberta::ElementInfo< dimension > &elementInfo
        = Grid::getRealImplementation( *it ).elementInfo();
      Alberta::ForLoop< MarkSubEntities, 2, dimension >::apply( hIndexSet, marker_, elementInfo );
    }
    up2Date_ = true;
  }


  template< int dim, int dimworld >
  inline void AlbertaMarkerVector< dim, dimworld >::print ( std::ostream &out ) const
  {
    for( int codim = 1; codim <= dimension; ++codim )
    {
      std::vector< int > &marker = marker_[ codim ];
      const int size = marker.size();
      if( size > 0)
      {
        out << std::endl;
        out << "Codimension " << codim << " (" << size << " entries)" << std::endl;
        for( int i = 0; i < size; ++i )
          out << "subentity " << i << " visited on Element " << marker[ i ] << std::endl;
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
                              const MarkerVector *vertexMark,
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
