// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_ITERATOR_IMP_CC
#define DUNE_ALUGRID_ITERATOR_IMP_CC

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include "alu3dinclude.hh"

#include "geometry.hh"
#include "entity.hh"
#include "grid.hh"
#include "faceutility.hh"

namespace Dune {

  /************************************************************************************
  ###
  #     #    #   #####  ######  #####    ####   ######   ####      #     #####
  #     ##   #     #    #       #    #  #       #       #    #     #       #
  #     # #  #     #    #####   #    #   ####   #####   #          #       #
  #     #  # #     #    #       #####        #  #       #          #       #
  #     #   ##     #    #       #   #   #    #  #       #    #     #       #
  ###    #    #     #    ######  #    #   ####   ######   ####      #       #
  ************************************************************************************/

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridIntersectionIterator<GridImp> ::
  ALU3dGridIntersectionIterator(const FactoryType& factory,
                                int wLevel) :
    connector_( factory.grid().conformingRefinement(), factory.grid().ghostCellsEnabled() ),
    geoProvider_(connector_),
    factory_( factory ),
    item_(0),
    ghost_(0),
    index_(-1)
  {}

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridIntersectionIterator<GridImp> ::
  ALU3dGridIntersectionIterator(const FactoryType& factory,
                                HElementType *el,
                                int wLevel,bool end) :
    connector_( factory.grid().conformingRefinement(), factory.grid().ghostCellsEnabled() ),
    geoProvider_(connector_),
    factory_( factory ),
    item_(0),
    ghost_(0),
    index_(-1)
  {
    if( ! end )
    {
      setFirstItem(*el,wLevel);
    }
    else
    {
      done();
    }
  }

  template<class GridImp>
  inline void
  ALU3dGridIntersectionIterator<GridImp> :: done ()
  {
    item_  = 0;
    ghost_ = 0;
    // index < 0 indicates end iterator
    index_ = -1;
  }

  template<class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp> ::
  setFirstItem (const HElementType & elem, int wLevel)
  {
    ghost_      = 0;
    item_       = static_cast<const IMPLElementType *> (&elem);

    // Get first face
    const GEOFaceType* firstFace = getFace(*item_, index_);

    const GEOFaceType* childFace = firstFace->down();
    if( childFace ) firstFace = childFace;

    // Store the face in the connector
    setNewFace(*firstFace);
  }

  template<class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp> ::
  setInteriorItem (const HElementType & elem, const BNDFaceType& ghost, int wLevel)
  {
    // get correct face number
    index_ = ElementTopo::alu2duneFace( ghost.getGhost().second );

    // store ghost for method inside
    ghost_   = &ghost;

    // Get first face
    const GEOFaceType* firstFace = getFace( ghost, index_ );
    item_   = static_cast<const IMPLElementType *> (&elem);

    const GEOFaceType* childFace = firstFace->down();
    if( childFace ) firstFace = childFace;

    // Store the face in the connector
    setGhostFace(*firstFace);
  }

  template<class GridImp>
  template <class EntityType>
  inline void ALU3dGridIntersectionIterator<GridImp> ::
  first (const EntityType & en, int wLevel)
  {
    if( ! en.isLeaf() )
    {
      done();
      return ;
    }

    innerLevel_ = en.level();
    index_  = 0;

    if( en.isGhost() )
    {
      setInteriorItem(en.getItem(), en.getGhost(), wLevel);
    }
    else
    {
      assert( numFaces == en.getItem().nFaces() );
      setFirstItem(en.getItem(), wLevel);
    }
  }

  // copy constructor
  template<class GridImp>
  inline ALU3dGridIntersectionIterator<GridImp> ::
  ALU3dGridIntersectionIterator(const ALU3dGridIntersectionIterator<GridImp> & org) :
    connector_(org.connector_),
    geoProvider_(connector_),
    factory_( org.factory_ ),
    item_(org.item_),
    ghost_(org.ghost_)
  {
    if(org.item_)
    { // else it's a end iterator
      item_        = org.item_;
      innerLevel_  = org.innerLevel_;
      index_       = org.index_;
    }
    else
    {
      done();
    }
  }

  // copy constructor
  template<class GridImp>
  inline void
  ALU3dGridIntersectionIterator<GridImp> ::
  assign(const ALU3dGridIntersectionIterator<GridImp> & org)
  {
    if(org.item_)
    {
      // else it's a end iterator
      item_       = org.item_;
      ghost_      = org.ghost_;
      innerLevel_ = org.innerLevel_;
      index_      = org.index_;
      connector_.updateFaceInfo(org.connector_.face(),innerLevel_,
                                item_->twist(ElementTopo::dune2aluFace(index_)));
      geoProvider_.resetFaceGeom();
    }
    else {
      done();
    }
    assert( equals(org) );
  }

  // check whether entities are the same or whether iterator is done
  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp> ::
  equals (const ALU3dGridIntersectionIterator<GridImp> & i ) const
  {
    // this method is only to check equality of real iterators and end
    // iterators
    return ((item_  == i.item_) &&
            (index_ == i.index_ )
            );
  }

  template<class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp> :: increment ()
  {
    // leaf increment
    assert(item_);

    const GEOFaceType * nextFace = 0;

    // When neighbour element is refined, try to get the next child on the face
    if (connector_.conformanceState() == FaceInfoType::REFINED_OUTER)
    {
      nextFace = connector_.face().next();

      // There was a next child face...
      if (nextFace)
      {
        if( ImplTraits :: isGhost( ghost_ ) )
        {
          setGhostFace( *nextFace );
        }
        else
        {
          setNewFace(*nextFace);
        }
        return; // we found what we were looking for...
      }
    } // end if

    // Next face number of starting element
    ++index_;

    // When the face number is larger than the number of faces an element
    // can have, we've reached the end...
    // for ghost elements here is finito
    if (index_ >= numFaces || ghost_ )
    {
      done();
      return;
    }

    // ... else we can take the next face
    nextFace = getFace(connector_.innerEntity(), index_);
    assert(nextFace);

    // Check whether we need to go down first
    //if (nextFace has children which need to be visited)
    const GEOFaceType * childFace = nextFace->down();
    if( childFace ) nextFace = childFace;

    assert(nextFace);
    setNewFace(*nextFace);
    return;
  }


  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::EntityPointer
  ALU3dGridIntersectionIterator<GridImp>::outside () const
  {
    assert( neighbor() );
    // make sure that outside is not called for an end iterator

    if( connector_.ghostBoundary() )
    {
      // create entity pointer with ghost boundary face
      return EntityPointer(factory_, connector_.boundaryFace() );
    }

    assert( &connector_.outerEntity() );
    return EntityPointer(factory_, connector_.outerEntity() );
  }

  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::EntityPointer
  ALU3dGridIntersectionIterator<GridImp>::inside () const
  {
    if( ImplTraits :: isGhost( ghost_ ) )
    {
      return EntityPointer(factory_, *ghost_ );
    }
    else
    {
      // make sure that inside is not called for an end iterator
      return EntityPointer(factory_, connector_.innerEntity() );
    }
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp> :: boundary () const
  {
    return connector_.boundary();
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp> :: neighbor () const
  {
    return connector_.neighbor();
  }

  template<class GridImp>
  inline int
  ALU3dGridIntersectionIterator< GridImp >::indexInInside () const
  {
    assert(ElementTopo::dune2aluFace(index_) == connector_.innerALUFaceIndex());
    return index_;
  }

  template< class GridImp >
  inline typename ALU3dGridIntersectionIterator< GridImp >::LocalGeometry
  ALU3dGridIntersectionIterator< GridImp >::geometryInInside () const
  {
    buildLocalGeometries();
    return LocalGeometry( intersectionSelfLocal_ );
  }


  template< class GridImp >
  inline int
  ALU3dGridIntersectionIterator< GridImp >::indexInOutside () const
  {
    return ElementTopo::alu2duneFace( connector_.outerALUFaceIndex() );
  }

  template< class GridImp >
  inline int
  ALU3dGridIntersectionIterator< GridImp >::
  twistInInside () const
  {
    return connector_.duneTwist( indexInInside(), connector_.innerTwist() );
  }

  template< class GridImp >
  inline int
  ALU3dGridIntersectionIterator< GridImp >::
  twistInOutside () const
  {
    return connector_.duneTwist( indexInOutside(), connector_.outerTwist() );
  }

  template< class GridImp >
  inline typename ALU3dGridIntersectionIterator< GridImp >::LocalGeometry
  ALU3dGridIntersectionIterator< GridImp >::geometryInOutside () const
  {
    assert(neighbor());
    buildLocalGeometries();
    return LocalGeometry( intersectionNeighborLocal_ );
  }

  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::NormalType &
  ALU3dGridIntersectionIterator<GridImp>::
  integrationOuterNormal(const FieldVector<alu3d_ctype, dim-1>& local) const
  {
    return this->outerNormal(local);
  }

  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::NormalType &
  ALU3dGridIntersectionIterator<GridImp>::
  outerNormal(const FieldVector<alu3d_ctype, dim-1>& local) const
  {
    assert(item_ != 0);
    return geoProvider_.outerNormal(local);
  }

  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::NormalType &
  ALU3dGridIntersectionIterator<GridImp>::
  unitOuterNormal(const FieldVector<alu3d_ctype, dim-1>& local) const
  {
    unitOuterNormal_ = this->outerNormal(local);
    unitOuterNormal_ *= (1.0/unitOuterNormal_.two_norm());
    return unitOuterNormal_;
  }

  template< class GridImp >
  inline typename ALU3dGridIntersectionIterator< GridImp >::Geometry
  ALU3dGridIntersectionIterator< GridImp >::geometry () const
  {
    geoProvider_.buildGlobalGeom( intersectionGlobal_ );
    return Geometry( intersectionGlobal_ );
  }

  template<class GridImp>
  inline GeometryType
  ALU3dGridIntersectionIterator<GridImp>::
  type () const
  {
    return GeometryType(
             GridImp::elementType == tetra ?
             GenericGeometry :: SimplexTopology< dim-1 > :: type :: id :
             GenericGeometry :: CubeTopology   < dim-1 > :: type :: id,
             dim-1 );
  }

  template<class GridImp>
  inline int
  ALU3dGridIntersectionIterator<GridImp>::boundaryId () const
  {
    assert( item_ );
    return ( boundary() ) ? connector_.boundaryId() : 0;
  }

  template<class GridImp>
  inline size_t
  ALU3dGridIntersectionIterator<GridImp>::boundarySegmentIndex() const
  {
    assert( item_ );
    assert( boundary() );
    return connector_.segmentIndex();
  }

  template< class GridImp >
  inline void ALU3dGridIntersectionIterator< GridImp >::buildLocalGeometries() const
  {
    intersectionSelfLocal_.buildGeom( geoProvider_.intersectionSelfLocal() );
    if ( !connector_.outerBoundary() )
      intersectionNeighborLocal_.buildGeom( geoProvider_.intersectionNeighborLocal() );
  }

  template <class GridImp>
  inline const typename ALU3dImplTraits< tetra, typename GridImp::MPICommunicatorType >::GEOFaceType *
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOTriangleBndType & bnd, int index) const
  {
    return bnd.myhface3(0);
  }

  template <class GridImp>
  inline const typename ALU3dImplTraits< hexa, typename GridImp::MPICommunicatorType >::GEOFaceType *
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOQuadBndType & bnd, int index) const
  {
    return bnd.myhface4(0);
  }

  template <class GridImp>
  inline const typename ALU3dImplTraits< tetra, typename GridImp::MPICommunicatorType >::GEOFaceType *
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOTetraElementType & elem, int index) const {
    assert(index >= 0 && index < numFaces);
    return elem.myhface3(ElementTopo::dune2aluFace(index));
  }

  template <class GridImp>
  inline const typename ALU3dImplTraits< hexa, typename GridImp::MPICommunicatorType >::GEOFaceType *
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOHexaElementType & elem, int index) const {
    assert(index >= 0 && index < numFaces);
    return elem.myhface4(ElementTopo::dune2aluFace(index));
  }

  template <class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp>::
  setNewFace(const GEOFaceType& newFace)
  {
    assert( ! ghost_ );
    assert( innerLevel_ == item_->level() );
    connector_.updateFaceInfo(newFace,innerLevel_,
                              item_->twist(ElementTopo::dune2aluFace(index_)) );
    geoProvider_.resetFaceGeom();
  }

  template <class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp>::
  setGhostFace(const GEOFaceType& newFace)
  {
    assert( ghost_ );
    assert( innerLevel_ == ghost_->level() );
    connector_.updateFaceInfo(newFace,innerLevel_, ghost_->twist(0) );
    geoProvider_.resetFaceGeom();
  }

  template <class GridImp>
  inline int
  ALU3dGridIntersectionIterator<GridImp>::
  level() const {
    assert( item_ && (innerLevel_ == item_->level()) );
    return innerLevel_;
  }

  /************************************************************************************
  ###
  #     #    #   #####  ######  #####    ####   ######   ####      #     #####
  #     ##   #     #    #       #    #  #       #       #    #     #       #
  #     # #  #     #    #####   #    #   ####   #####   #          #       #
  #     #  # #     #    #       #####        #  #       #          #       #
  #     #   ##     #    #       #   #   #    #  #       #    #     #       #
  ###    #    #     #    ######  #    #   ####   ######   ####      #       #
  ************************************************************************************/

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridLevelIntersectionIterator<GridImp> ::
  ALU3dGridLevelIntersectionIterator(const FactoryType& factory,
                                     int wLevel)
    : ALU3dGridIntersectionIterator<GridImp>(factory,wLevel)
      , levelNeighbor_(false)
      , isLeafItem_(false)
  {}

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridLevelIntersectionIterator<GridImp> ::
  ALU3dGridLevelIntersectionIterator(const FactoryType& factory,
                                     HElementType *el,
                                     int wLevel,bool end)
    : ALU3dGridIntersectionIterator<GridImp>(factory,el,wLevel,end)
      , levelNeighbor_(false)
      , isLeafItem_(false)
  {}

  template<class GridImp>
  template <class EntityType>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> ::
  first (const EntityType & en, int wLevel)
  {
    // if given Entity is not leaf, we create an end iterator
    index_  = 0;
    isLeafItem_   = en.isLeaf();

    if( en.isGhost() )
    {
      setInteriorItem(en.getItem(), en.getGhost(), wLevel);
    }
    else
    {
      assert( numFaces == en.getItem().nFaces() );
      setFirstItem(en.getItem(), wLevel);
    }
  }

  template<class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> ::
  setFirstItem (const HElementType & elem, int wLevel)
  {
    ghost_       = 0;
    item_        = static_cast<const IMPLElementType *> (&elem);
    this->innerLevel_  = wLevel;
    // Get first face
    const GEOFaceType* firstFace = getFace(*item_, index_);
    // Store the face in the connector
    setNewFace(*firstFace);
  }

  template<class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> ::
  setInteriorItem (const HElementType & elem, const BNDFaceType& ghost, int wLevel)
  {
    // store ghost for method inside
    ghost_   = &ghost;
    item_   = static_cast<const IMPLElementType *> (&elem);
    // get correct face number
    index_ = ElementTopo::alu2duneFace( ghost.getGhost().second );

    innerLevel_  = wLevel;

    // Get first face
    const GEOFaceType* firstFace = getFace( ghost, index_ );

    // Store the face in the connector
    setNewFace(*firstFace);
  }

  // copy constructor
  template<class GridImp>
  inline ALU3dGridLevelIntersectionIterator<GridImp> ::
  ALU3dGridLevelIntersectionIterator(const ThisType & org)
    : ALU3dGridIntersectionIterator<GridImp>(org)
      , levelNeighbor_(org.levelNeighbor_)
      , isLeafItem_(org.isLeafItem_)
  {}

  // copy constructor
  template<class GridImp>
  inline void
  ALU3dGridLevelIntersectionIterator<GridImp> ::
  assign(const ALU3dGridLevelIntersectionIterator<GridImp> & org)
  {
    ALU3dGridIntersectionIterator<GridImp>::assign(org);
    levelNeighbor_ = org.levelNeighbor_;
    isLeafItem_    = org.isLeafItem_;
  }

  template<class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> :: increment ()
  {
    // level increment
    assert( item_ );

    // Next face number of starting element
    ++index_;

    // When the face number is larger than the number of faces an element
    // can have, we've reached the end...
    if ( index_ >= numFaces || ImplTraits::isGhost( ghost_ ) )
    {
      done();
      return;
    }

    // ... else we can take the next face
    const GEOFaceType * nextFace = getFace(connector_.innerEntity(), index_);
    assert(nextFace);

    setNewFace(*nextFace);
    return;
  }

  template<class GridImp>
  inline bool ALU3dGridLevelIntersectionIterator<GridImp>::neighbor () const
  {
    return levelNeighbor_ && (BaseType :: neighbor());
  }

  template <class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp>::
  setNewFace(const GEOFaceType& newFace)
  {
    assert( item_->level() == innerLevel_ );
    levelNeighbor_ = (newFace.level() == innerLevel_);
    connector_.updateFaceInfo(newFace, innerLevel_,
                              ( ImplTraits::isGhost( ghost_ ) ) ?
                              ghost_->twist(0) :
                              item_->twist(ElementTopo::dune2aluFace( index_ )));
    geoProvider_.resetFaceGeom();

    // check again level neighbor because outer element might be coarser then
    // this element
    if( isLeafItem_ )
    {
      if( connector_.ghostBoundary() )
      {
        const BNDFaceType & ghost = connector_.boundaryFace();
        // if nonconformity occurs then no level neighbor
        levelNeighbor_ = (innerLevel_ == ghost.ghostLevel() );
      }
      else if ( ! connector_.outerBoundary() )
      {
        levelNeighbor_ = (connector_.outerEntity().level() == innerLevel_);
      }
    }
  }

} // end namespace Dune

#if COMPILE_ALUGRID_INLINE
  #include "iterator.cc"
#endif

#endif // DUNE_ALUGRID_ITERATOR_IMP_CC
