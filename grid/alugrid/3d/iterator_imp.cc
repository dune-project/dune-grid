// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "geometry.hh"
#include "entity.hh"
#include "grid.hh"
#include "faceutility.hh"
#include "alu3dutility.hh"

#include <dune/common/misc.hh>

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
  ALU3dGridIntersectionIterator(const GridImp & grid,
                                int wLevel) :
    geoProvider_(connector_),
    grid_(grid),
    item_(0),
    index_(0),
    intersectionGlobal_(GeometryImp()),
    intersectionGlobalImp_(grid_.getRealImplementation(intersectionGlobal_)),
    intersectionSelfLocal_(GeometryImp()),
    intersectionSelfLocalImp_(grid_.getRealImplementation(intersectionSelfLocal_)),
    intersectionNeighborLocal_(GeometryImp()),
    intersectionNeighborLocalImp_(grid_.getRealImplementation(intersectionNeighborLocal_)),
    done_(true)
  {}

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridIntersectionIterator<GridImp> ::
  ALU3dGridIntersectionIterator(const GridImp & grid,
                                HElementType *el,
                                int wLevel,bool end) :
    connector_(),
    geoProvider_(connector_),
    grid_(grid),
    item_(0),
    index_(0),
    intersectionGlobal_(GeometryImp()),
    intersectionGlobalImp_(grid_.getRealImplementation(intersectionGlobal_)),
    intersectionSelfLocal_(GeometryImp()),
    intersectionSelfLocalImp_(grid_.getRealImplementation(intersectionSelfLocal_)),
    intersectionNeighborLocal_(GeometryImp()),
    intersectionNeighborLocalImp_(grid_.getRealImplementation(intersectionNeighborLocal_)),
    done_(end)
  {
    if (!end)
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
    done_ = true;
    item_ = 0;
  }

  template<class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp> ::
  setFirstItem (const HElementType & elem, int wLevel)
  {
    item_       = static_cast<const IMPLElementType *> (&elem);

    // Get first face
    const GEOFaceType* firstFace = getFace(*item_, index_);

    const GEOFaceType * childFace = firstFace->down();
    if( childFace ) firstFace = childFace;

    // Store the face in the connector
    setNewFace(*firstFace);
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

    done_   = false;
    assert( numFaces == en.getItem().nFaces() );
    innerLevel_ = en.level();
    index_  = 0;
    setFirstItem(en.getItem(),wLevel);
  }

  // copy constructor
  template<class GridImp>
  inline ALU3dGridIntersectionIterator<GridImp> ::
  ALU3dGridIntersectionIterator(const ALU3dGridIntersectionIterator<GridImp> & org) :
    connector_(org.connector_),
    geoProvider_(connector_),
    grid_(org.grid_),
    item_(org.item_),
    intersectionGlobal_(GeometryImp()),
    intersectionGlobalImp_(grid_.getRealImplementation(intersectionGlobal_)),
    intersectionSelfLocal_(GeometryImp()),
    intersectionSelfLocalImp_(grid_.getRealImplementation(intersectionSelfLocal_)),
    intersectionNeighborLocal_(GeometryImp()),
    intersectionNeighborLocalImp_(grid_.getRealImplementation(intersectionNeighborLocal_)),
    done_(org.done_)
  {
    if(org.item_) { // else it's a end iterator
      item_        = org.item_;
      innerLevel_  = org.innerLevel_;
      index_       = org.index_;
    } else {
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
      item_      = org.item_;
      innerLevel_ = org.innerLevel_;
      index_     = org.index_;
      connector_.updateFaceInfo(org.connector_.face(),innerLevel_,
                                item_->twist(ElementTopo::dune2aluFace(index_)));
      geoProvider_.resetFaceGeom();
    }
    else {
      done();
    }
  }

  // check whether entities are the same or whether iterator is done
  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp> ::
  equals (const ALU3dGridIntersectionIterator<GridImp> & i ) const
  {
    // this method is only to check equality of real iterators and end
    // iterators
    return ((item_ == i.item_) &&
            (done_ == i.done_)
            );
  }

  template<class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp> :: increment ()
  {
    // leaf increment
    assert(item_);

    const GEOFaceType * nextFace = 0;

    // When neighbour element is refined, try to get the next child on the face
    if (connector_.conformanceState() == FaceInfoType::REFINED_OUTER) {
      nextFace = connector_.face().next();

      // There was a next child face...
      if (nextFace) {
        setNewFace(*nextFace);
        return; // we found what we were looking for...
      }
    } // end if

    // Next face number of starting element
    ++index_;

    // When the face number is larger than the number of faces an element
    // can have, we've reached the end...
    if (index_ >= numFaces) {
      this->done();
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
    assert( this->neighbor() );
    // make sure that outside is not called for an end iterator
#if ALU3DGRID_PARALLEL
    if(connector_.ghostBoundary())
    {
      // create entity pointer with ghost boundary face
      return EntityPointer(this->grid_, connector_.boundaryFace() );
    }
#endif
    assert( &connector_.outerEntity() );
    return EntityPointer(this->grid_, connector_.outerEntity() );
  }

  template<class GridImp>
  inline typename ALU3dGridIntersectionIterator<GridImp>::EntityPointer
  ALU3dGridIntersectionIterator<GridImp>::inside () const {
    // make sure that inside is not called for an end iterator
    //assert( !done_ );
    return EntityPointer(this->grid_, connector_.innerEntity() );
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp> :: boundary () const
  {
    return (connector_.outerBoundary());
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp>::neighbor () const
  {
    return ! boundary();
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp>::levelNeighbor () const
  {
    return false;
  }

  template<class GridImp>
  inline bool ALU3dGridIntersectionIterator<GridImp>::leafNeighbor () const
  {
    return neighbor();
  }


  template<>
  inline int
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, tetra > >::indexInInside () const
  {
    assert(ElementTopo::dune2aluFace(index_) == connector_.innerALUFaceIndex());
    return 3 - index_;
  }

  template<>
  inline int
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, hexa > >::indexInInside () const
  {
    assert(ElementTopo::dune2aluFace(index_) == connector_.innerALUFaceIndex());
    return index_;
  }


  template <class GridImp>
  inline const typename ALU3dGridIntersectionIterator<GridImp>::LocalGeometry &
  ALU3dGridIntersectionIterator< GridImp >::geometryInInside () const
  {
    buildLocalGeometries();
    return intersectionSelfLocal_;
  }


  template<>
  inline int
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, tetra > >::indexInOutside () const
  {
    return 3 - ElementTopo::alu2duneFace( connector_.outerALUFaceIndex() );
  }

  template<>
  inline int
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, hexa > >::indexInOutside () const
  {
    return ElementTopo::alu2duneFace( connector_.outerALUFaceIndex() );
  }


  template<>
  inline int ALU3dGridIntersectionIterator< const ALU3dGrid<3,3,hexa> >::twistInSelf () const
  {
    const int aluTwist = connector_.innerTwist();
    const int mappedZero =
      FaceTopo::twist(ElementTopo::dune2aluFaceVertex( indexInInside(), 0), aluTwist);

    return
      (ElementTopo::faceOrientation( indexInInside() ) * sign(aluTwist) < 0 ?
       mappedZero : -mappedZero-1);
  }

  template<>
  inline int ALU3dGridIntersectionIterator< const ALU3dGrid<3 ,3, tetra> >::twistInSelf () const
  {
    return connector_.innerTwist();
  }

  template<>
  inline int ALU3dGridIntersectionIterator< const ALU3dGrid<3,3,hexa> >::twistInNeighbor () const
  {
    const int aluTwist = connector_.outerTwist();
    const int mappedZero =
      FaceTopo::twist(ElementTopo::dune2aluFaceVertex( indexInOutside(), 0), aluTwist);

    return
      (ElementTopo::faceOrientation( indexInOutside() ) * sign(aluTwist) < 0 ?
       mappedZero : -mappedZero-1);
  }

  template<>
  inline int ALU3dGridIntersectionIterator< const ALU3dGrid<3 ,3, tetra> >::twistInNeighbor () const
  {
    return connector_.outerTwist();
  }

  template <class GridImp>
  inline const typename ALU3dGridIntersectionIterator<GridImp>::LocalGeometry &
  ALU3dGridIntersectionIterator< GridImp >::geometryInOutside () const
  {
    assert(neighbor());
    buildLocalGeometries();
    return intersectionNeighborLocal_;
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

  template<class GridImp>
  inline const typename ALU3dGridIntersectionIterator<GridImp>::Geometry &
  ALU3dGridIntersectionIterator< GridImp >::geometry () const
  {
    geoProvider_.buildGlobalGeom(intersectionGlobalImp_);
    return intersectionGlobal_;
  }


  template<>
  inline GeometryType
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, tetra > >::type () const
  {
    return GeometryType( GeometryType::simplex, dim-1 );
  }

  template<>
  inline GeometryType
  ALU3dGridIntersectionIterator< const ALU3dGrid< 3, 3, hexa > >::type () const
  {
    return GeometryType( GeometryType::cube, dim-1 );
  }


  template<class GridImp>
  inline int
  ALU3dGridIntersectionIterator<GridImp>::boundaryId () const
  {
    assert(item_);
    return (boundary() ? connector_.boundaryFace().bndtype() : 0);
  }

  template <class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp>::buildLocalGeometries() const
  {
    intersectionSelfLocalImp_.buildGeom(geoProvider_.intersectionSelfLocal());
    if (!connector_.outerBoundary()) {
      intersectionNeighborLocalImp_.buildGeom(geoProvider_.intersectionNeighborLocal());
    }
  }

  template <class GridImp>
  inline const ALU3dImplTraits<tetra>::GEOFaceType*
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOTetraElementType & elem, int index) const {
    assert(index >= 0 && index < numFaces);
    return elem.myhface3(ElementTopo::dune2aluFace(index));
  }

  template <class GridImp>
  inline const ALU3dImplTraits<hexa>::GEOFaceType*
  ALU3dGridIntersectionIterator<GridImp>::
  getFace(const GEOHexaElementType & elem, int index) const {
    assert(index >= 0 && index < numFaces);
    return elem.myhface4(ElementTopo::dune2aluFace(index));
  }

  template <class GridImp>
  inline void ALU3dGridIntersectionIterator<GridImp>::
  setNewFace(const GEOFaceType& newFace)
  {
    assert( innerLevel_ == item_->level() );
    connector_.updateFaceInfo(newFace,innerLevel_,
                              item_->twist(ElementTopo::dune2aluFace(index_)));
    geoProvider_.resetFaceGeom();
  }

  template <class GridImp>
  void ALU3dGridIntersectionIterator<GridImp>::outputElementInfo() const
  {
    std::cout << "Starting outputElementInfo\n";
    // output element corner coordinates
    std::cout << "Element corner coordinates" << std::endl;
    for (int i = 0; i < numVertices; ++i) {
      printToScreen(ElementTopo::alu2duneVertex(i), i,
                    convert2FV(item_->myvertex(i)->Point()));
    }

    // output midpoint of faces
    /* too specific for hexa
       std::cout << "Midpoint of faces" << std::endl;
       FieldVector<alu3d_ctype, 2> mid(0.5);
       for (int i = 0; i < numFaces; ++i) {
       NormalType c0 =
        convert2FV(item_->myvertex(i, FaceTopo::dune2aluVertex(0))->Point());
       NormalType c1 =
        convert2FV(item_->myvertex(i, FaceTopo::dune2aluVertex(1))->Point());
       NormalType c2 =
        convert2FV(item_->myvertex(i, FaceTopo::dune2aluVertex(2))->Point());
       NormalType c3 =
        convert2FV(item_->myvertex(i, FaceTopo::dune2aluVertex(3))->Point());

       BilinearSurfaceMapping biMap(c0, c1, c2, c3);
       NormalType result;
       biMap.map2world(mid, result);

       printToScreen(ElementTopo::alu2duneFace(i), i, result);
       } // end for

       // output element reference faces - as calculated in Dune
       std::cout << "Corner indices (from ReferenceTopologySet)" << std::endl;
       for (int i = 0; i < numFaces; ++i) {
       std::cout << "-- Face " << i << "--\n";
       const int* result;
       int n;
       ReferenceTopologySet::getSubEntities< 1, 3 >(hexahedron,
                                                   i,
                                                   result,
                                                   n);
       for (int j = 0; j < n; ++j) {
        printToScreen(result[j], ElementTopo::dune2aluVertex(result[j]), j);
       } // end for
       } // end for
     */
    // output element reference faces - compared with data from Alu
    std::cout << "Corner indices (from ALU)" << std::endl;
    for (int i = 0; i < numFaces; ++i) {
      std::cout << "-- Face " << i << "--\n";
      for (int j = 0; j < numVerticesPerFace; ++j) {
        int aluIdx = GEOElementType::prototype[i][j];
        printToScreen(ElementTopo::alu2duneVertex(aluIdx), aluIdx, j);
      }
    }

    std::cout << "Ending outputElementInfo\n" << std::endl;
  }

  template <class GridImp>
  void ALU3dGridIntersectionIterator<GridImp>::outputFaceInfo() const {
    std::cout << "Starting outputFaceInfo\n";
    // output face index (inner, outer)
    std::cout << "Inner twist" << std::endl;
    printToScreen(index_,
                  connector_.innerALUFaceIndex(),
                  connector_.innerTwist());

    assert(index_ == ElementTopo::alu2duneFace(connector_.innerALUFaceIndex()));

    std::cout << "Outer twist" << std::endl;
    printToScreen(-1,
                  connector_.outerALUFaceIndex(),
                  connector_.outerTwist());

    // output corner coordinates
    std::cout << "Face corner coordinates (ALU face)" << std::endl;
    const GEOFaceType& face = connector_.face();
    for (int i = 0; i < numVerticesPerFace; ++i) {
      printToScreen(FaceTopo::alu2duneVertex(i), i,
                    convert2FV(face.myvertex(i)->Point()));
    }

    std::cout << "Face corner coordinates (intersectionGlobal)" << std::endl;
    const Geometry& interGlobal = geometry();
    for (int i = 0; i < numVerticesPerFace; ++i) {
      printToScreen(i, FaceTopo::dune2aluVertex(i), interGlobal[i]);
    }

    std::cout << "Ending outputFaceInfo\n" << std::endl;
  }

  template <class GridImp>
  void ALU3dGridIntersectionIterator<GridImp>::
  printToScreen(int duneIdx, int aluIdx) const {
    printToScreen(duneIdx, aluIdx, "-");
  }

  template <class GridImp>
  template <typename T>
  void ALU3dGridIntersectionIterator<GridImp>::
  printToScreen(int duneIdx, int aluIdx, const T& info) const {
    std::cout << duneIdx << ", " << aluIdx << ": " << info << "\n";
  }

  template <class GridImp>
  typename ALU3dGridIntersectionIterator<GridImp>::NormalType
  ALU3dGridIntersectionIterator<GridImp>::
  convert2FV(const alu3d_ctype (&p)[3]) const {
    NormalType result;
    result[0] = p[0];
    result[1] = p[1];
    result[2] = p[2];
    return result;
  }

  template <class GridImp>
  inline int
  ALU3dGridIntersectionIterator<GridImp>::
  level() const {
    assert( item_ );
    return item_->level();
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
  ALU3dGridLevelIntersectionIterator(const GridImp & grid,
                                     int wLevel)
    : ALU3dGridIntersectionIterator<GridImp>(grid,wLevel)
      , levelNeighbor_(false)
      , isLeafItem_(false)
  {}

  // --IntersectionIterator
  template<class GridImp>
  inline ALU3dGridLevelIntersectionIterator<GridImp> ::
  ALU3dGridLevelIntersectionIterator(const GridImp & grid,
                                     HElementType *el,
                                     int wLevel,bool end)
    : ALU3dGridIntersectionIterator<GridImp>(grid,el,wLevel,end)
      , levelNeighbor_(false)
      , isLeafItem_(false)
  {}

  template<class GridImp>
  template <class EntityType>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> ::
  first (const EntityType & en, int wLevel)
  {
    // if given Entity is not leaf, we create an end iterator
    this->done_   = false;
    assert( numFaces == en.getItem().nFaces() );
    this->index_  = 0;
    isLeafItem_ = en.isLeaf();
    setFirstItem(en.getItem(),wLevel);
  }

  template<class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp> ::
  setFirstItem (const HElementType & elem, int wLevel)
  {
    this->item_        = static_cast<const IMPLElementType *> (&elem);
    this->innerLevel_  = wLevel;
    // Get first face
    const GEOFaceType* firstFace = getFace(*this->item_, this->index_);
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
    assert(this->item_);

    // Next face number of starting element
    ++this->index_;

    // When the face number is larger than the number of faces an element
    // can have, we've reached the end...
    if (this->index_ >= numFaces) {
      this->done();
      return;
    }

    // ... else we can take the next face
    const GEOFaceType * nextFace = this->getFace(this->connector_.innerEntity(), this->index_);
    assert(nextFace);

    setNewFace(*nextFace);
    return;
  }

  template<class GridImp>
  inline bool ALU3dGridLevelIntersectionIterator<GridImp>::neighbor () const
  {
    return levelNeighbor();
  }

  template<class GridImp>
  inline bool ALU3dGridLevelIntersectionIterator<GridImp>::levelNeighbor () const
  {
    return levelNeighbor_ && (! this->boundary());
  }

  template<class GridImp>
  inline bool ALU3dGridLevelIntersectionIterator<GridImp>::leafNeighbor () const
  {
    return false;
  }


  template <class GridImp>
  inline void ALU3dGridLevelIntersectionIterator<GridImp>::
  setNewFace(const GEOFaceType& newFace)
  {
    assert( this->item_->level() == this->innerLevel_ );
    //const int itemLevel = wLevel; //this->item_->level();
    levelNeighbor_ = (newFace.level() == this->innerLevel_);
    this->connector_.updateFaceInfo(newFace,this->innerLevel_,
                                    this->item_->twist(ElementTopo::dune2aluFace(this->index_)));
    this->geoProvider_.resetFaceGeom();

    // check again level neighbor because outer element might be coarser then
    // this element
    if( isLeafItem_ )
    {
      if( this->connector_.ghostBoundary() )
      {
        const BNDFaceType & ghost = this->connector_.boundaryFace();
        // if nonconformity occurs then no level neighbor
        levelNeighbor_ = (this->innerLevel_ == ghost.ghostLevel() );
      }
      else if ( ! this->connector_.outerBoundary() )
      {
        levelNeighbor_ = (this->connector_.outerEntity().level() == this->innerLevel_);
      }
    }
  }
  /*************************************************************************
  #       ######  #    #  ######  #          #     #####  ######  #####
  #       #       #    #  #       #          #       #    #       #    #
  #       #####   #    #  #####   #          #       #    #####   #    #
  #       #       #    #  #       #          #       #    #       #####
  #       #        #  #   #       #          #       #    #       #   #
  ######  ######    ##    ######  ######     #       #    ######  #    #
  *************************************************************************/
  //--LevelIterator
  // Constructor for begin iterator
  template<int codim, PartitionIteratorType pitype, class GridImp >
  inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const GridImp & grid, int level, bool )
    : ALU3dGridEntityPointer<codim,GridImp> (grid,level)
      , level_(level)
      , iter_ (0)
  {
    iter_  = new IteratorType ( this->grid_ , level_, grid.nlinks() );
    assert( iter_ );
    this->firstItem(this->grid_,*this,level_);
  }

  // Constructor for end iterator
  template<int codim, PartitionIteratorType pitype, class GridImp >
  inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const GridImp & grid, int level)
    : ALU3dGridEntityPointer<codim,GridImp> (grid ,level)
      , level_(level)
      , iter_ (0)
  {
    this->done();
  }

  template<int codim, PartitionIteratorType pitype, class GridImp >
  inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const ALU3dGridLevelIterator<codim,pitype,GridImp> & org )
    : ALU3dGridEntityPointer<codim,GridImp> ( org.grid_ , org.level_ )
      , level_( org.level_ )
      , iter_(0)
  {
    assign(org);
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  ~ALU3dGridLevelIterator ()
  {
    removeIter();
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline void ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  removeIter ()
  {
    this->done();
    if(iter_)
    {
      delete iter_;
      iter_ = 0;
    }
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline void ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  assign(const ThisType & org)
  {
    assert( iter_ == 0 );
    ALU3dGridEntityPointer <codim,GridImp> :: clone (org);
    level_ = org.level_;
    if( org.iter_ )
    {
      iter_ = new IteratorType ( *(org.iter_) );
      assert( iter_ );
      if(!(iter_->done()))
      {
        this->setItem(this->grid_, *this, *iter_, level_ );
        assert( this->equals(org) );
      }
    }
    else
    {
      this->done();
    }
  }
  template<int codim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLevelIterator<codim, pitype, GridImp> &
  ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    removeIter();
    assign(org);
    return *this;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp >
  inline void ALU3dGridLevelIterator<codim,pitype,GridImp> :: increment ()
  {
    this->incrementIterator(this->grid_,*this,level_);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline typename ALU3dGridLevelIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLevelIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLevelIterator<cdim, pitype, GridImp> endIterator (this->grid_,level_);
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }

  //*******************************************************************
  //
  //  LEAFITERATOR
  //
  //--LeafIterator
  //*******************************************************************
  // constructor for end iterators
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator( const GridImp &grid, int level )
    : ALU3dGridEntityPointer <cdim,GridImp> ( grid, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    this->done();
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const GridImp &grid, int level ,
                        bool isBegin)
    : ALU3dGridEntityPointer <cdim,GridImp> ( grid, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    // create interior iterator
    iter_ = new IteratorType ( this->grid_ , level , grid.nlinks() );
    assert( iter_ );
    // -1 to identify as leaf iterator
    this->firstItem(this->grid_,*this,-1);
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const ThisType & org)
    : ALU3dGridEntityPointer <cdim,GridImp> ( org.grid_, org.walkLevel_ )
      , iter_(0)
      , walkLevel_(org.walkLevel_)
  {
    // assign iterator without cloning entity pointer again
    assign(org);
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ~ALU3dGridLeafIterator()
  {
    removeIter();
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  removeIter ()
  {
    this->done();
    if(iter_)
    {
      delete iter_;
      iter_ = 0;
    }
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU3dGridLeafIterator<cdim, pitype, GridImp> &
  ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    removeIter();
    assign(org);
    return *this;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  assign (const ThisType & org)
  {
    assert( iter_ == 0 );
    ALU3dGridEntityPointer <cdim,GridImp> :: clone (org);

    if( org.iter_ )
    {
      assert( !org.iter_->done() );
      iter_ = new IteratorType ( *(org.iter_) );
      assert( iter_ );

      if( !(iter_->done() ))
      {
        assert( !iter_->done());
        assert( !org.iter_->done() );
        // -1 to identify leaf iterator
        this->setItem(this->grid_,*this, *iter_,-1);
        assert( this->equals(org) );
      }
    }
    else
    {
      this->done();
    }

    walkLevel_ = org.walkLevel_;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> :: increment ()
  {
    // -1 to identify leaf iterator
    this->incrementIterator(this->grid_,*this,-1);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline typename ALU3dGridLeafIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLeafIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLeafIterator<cdim, pitype, GridImp> endIterator (this->grid_, this->grid_.maxLevel());
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }


  /************************************************************************************
  #     #
  #     #     #    ######  #####      #     #####  ######  #####
  #     #     #    #       #    #     #       #    #       #    #
  #######     #    #####   #    #     #       #    #####   #    #
  #     #     #    #       #####      #       #    #       #####
  #     #     #    #       #   #      #       #    #       #   #
  #     #     #    ######  #    #     #       #    ######  #    #
  ************************************************************************************/
  // --HierarchicIterator
  template <class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const GridImp & grid ,
                              const HElementType & elem, int maxlevel ,bool end)
    : ALU3dGridEntityPointer<0,GridImp> ( grid, maxlevel )
      , elem_(&elem)
      , maxlevel_(maxlevel)
  {
    if (!end)
    {
      HElementType * item =
        const_cast<HElementType *> (elem.down());
      if(item)
      {
        // we have children and they lie in the disired level range
        if(item->level() <= maxlevel_)
        {
          this->updateEntityPointer( item );
        }
        else
        { // otherwise do nothing
          this->done();
        }
      }
      else
      {
        this->done();
      }
    }
  }

  template <class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const ThisType& org)
    : ALU3dGridEntityPointer<0,GridImp> ( org.grid_, org.maxlevel_ )
  {
    assign( org );
  }

  template <class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp> &
  ALU3dGridHierarchicIterator<GridImp> ::
  operator = (const ThisType& org)
  {
    assign( org );
    return *this;
  }

  template <class GridImp>
  inline void
  ALU3dGridHierarchicIterator<GridImp> ::
  assign(const ThisType& org)
  {
    elem_     = org.elem_;
    maxlevel_ = org.maxlevel_;

    // this method will free entity
    ALU3dGridEntityPointer<0,GridImp> :: clone(org);
  }

  template <class GridImp>
  inline typename ALU3dGridHierarchicIterator<GridImp>::HElementType*
  ALU3dGridHierarchicIterator<GridImp>::
  goNextElement(HElementType * oldelem )
  {
    // strategy is:
    // - go down as far as possible and then over all children
    // - then go to father and next and down again

    HElementType * nextelem = oldelem->down();
    if(nextelem)
    {
      if(nextelem->level() <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->next();
    if(nextelem)
    {
      if(nextelem->level() <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->up();
    if(nextelem == elem_) return 0;

    while( !nextelem->next() )
    {
      nextelem = nextelem->up();
      if(nextelem == elem_) return 0;
    }

    if(nextelem) nextelem = nextelem->next();

    return nextelem;
  }

  template <class GridImp>
  inline void ALU3dGridHierarchicIterator<GridImp> :: increment ()
  {
    assert(this->item_ != 0);

    HElementType * nextItem = goNextElement( this->item_ );
    if(!nextItem)
    {
      this->done();
      return ;
    }

    this->updateEntityPointer(nextItem);
    return ;
  }

  template <class GridImp>
  inline typename ALU3dGridHierarchicIterator<GridImp> :: Entity &
  ALU3dGridHierarchicIterator<GridImp> :: dereference () const
  {
    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }




} // end namespace Dune
