// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/common/exceptions.hh>
#include <dune/grid/common/referenceelements.hh>

#include "geometry.hh"
#include "grid.hh"

namespace Dune {

  // --Entity
  template <int cd, int dim, class GridImp>
  inline ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const GridImp  &grid, int level) :
    grid_(grid),
    level_(0),
    gIndex_(-1),
    twist_(0),
    face_(-1),
    item_(0),
    geo_( GeometryImp() ),
    geoImp_(grid.getRealImplementation(geo_)),
    builtgeometry_(false),
    partitionType_(InteriorEntity)
  {}

  // --Entity
  template <int cd, int dim, class GridImp>
  inline ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<cd,dim,GridImp> & org) :
    grid_(org.grid_),
    level_(org.level_),
    gIndex_(org.gIndex_),
    twist_(org.twist_),
    face_(org.face_),
    item_(org.item_),
    geo_( GeometryImp() ),
    geoImp_(grid_.getRealImplementation(geo_)),
    builtgeometry_(false),
    partitionType_(org.partitionType_)
  {}

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  reset( int l )
  {
    item_  = 0;
    level_ = l;
    twist_ = 0;
    face_  = -1;
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  removeElement()
  {
    item_ = 0;
  }

  template<int cd, int dim, class GridImp>
  inline bool ALU3dGridEntity<cd,dim,GridImp> ::
  equals(const ALU3dGridEntity<cd,dim,GridImp> & org) const
  {
    return (item_ == org.item_);
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setEntity(const ALU3dGridEntity<cd,dim,GridImp> & org)
  {
    item_   = org.item_;
    gIndex_ = org.gIndex_;
    twist_  = org.twist_;
    level_  = org.level_;
    face_   = org.face_;
    builtgeometry_= false;
    partitionType_ = org.partitionType_;
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const ElementType & item)
  {
    setElement(item,GetLevel<GridImp,cd>::getLevel(grid_,item));
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const ElementType & item, const int level, int twist , int face )
  {
    item_   = static_cast<const IMPLElementType *> (&item);
    gIndex_ = (*item_).getIndex();
    twist_  = twist;
    level_  = level;
    face_   = face;
    builtgeometry_=false;
    partitionType_ = this->convertBndId( *item_ );
  }

  template<>
  inline void ALU3dGridEntity<3,3,const ALU3dGrid<3,3,hexa> > ::
  setElement(const ALU3DSPACE HElementType &el, const ALU3DSPACE VertexType &vx)
  {
    item_   = static_cast<const IMPLElementType *> (&vx);
    gIndex_ = (*item_).getIndex();
    level_  = (*item_).level();
    builtgeometry_=false;
    partitionType_ = this->convertBndId( *item_ );
  }

  template<>
  inline void ALU3dGridEntity<3,3,const ALU3dGrid<3,3,tetra> > ::
  setElement(const ALU3DSPACE HElementType &el, const ALU3DSPACE VertexType &vx)
  {
    // * what the heck does this static_cast do!?
    item_   = static_cast<const IMPLElementType *> (&vx);
    gIndex_ = (*item_).getIndex();
    level_  = (*item_).level();
    builtgeometry_=false;
    partitionType_ = this->convertBndId( *item_ );
  }

  template<int cd, int dim, class GridImp>
  inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setGhost(const ALU3DSPACE HBndSegType &ghost)
  {
    // this method only exists, that we don't have to pecialise the
    // Iterators for each codim, this method should not be called otherwise
    // error
    DUNE_THROW(GridError,"This method should not be called!");
  }

  template<int cd, int dim, class GridImp>
  inline int ALU3dGridEntity<cd,dim,GridImp> :: getIndex () const
  {
    return gIndex_;
  }

  template<int cd, int dim, class GridImp>
  inline int ALU3dGridEntity<cd,dim,GridImp> :: level () const
  {
    return level_;
  }

  template<int cd, int dim, class GridImp>
  inline PartitionType ALU3dGridEntity<cd,dim,GridImp> ::
  partitionType () const
  {
    return partitionType_;
  }

  template<int cd, int dim, class GridImp>
  inline const typename ALU3dGridEntity<cd,dim,GridImp>::Geometry &
  ALU3dGridEntity<cd,dim,GridImp>:: geometry() const
  {
    //assert( (cd == 1) ? (face_ >= 0) : 1 );
    if(!builtgeometry_) builtgeometry_ = geoImp_.buildGeom(*item_, twist_, face_ );
    return geo_;
  }

  template<int cd, int dim, class GridImp>
  inline GeometryType
  ALU3dGridEntity<cd,dim,GridImp>:: type () const
  {
    return geo_.type();
  }

  template<int dim, class GridImp>
  inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const GridImp  &grid, int wLevel)
    : grid_(grid)
      , item_( 0 )
      , ghost_( 0 )
      , geo_(GeometryImp())
      , geoImp_ (grid_.getRealImplementation(geo_))
      , builtgeometry_(false)
      , walkLevel_ (wLevel)
      , level_(-1)
      , geoInFather_ (GeometryImp())
      , geoInFatherImp_(grid_.getRealImplementation(geoInFather_))
      , isLeaf_ (false)
      , refElem_(grid_.referenceElement())
  {  }

  template<int dim, class GridImp>
  inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<0,dim,GridImp> & org)
    : grid_(org.grid_)
      , item_(org.item_)
      , ghost_( org.ghost_ )
      , geo_(GeometryImp())
      , geoImp_ (grid_.getRealImplementation(geo_))
      , builtgeometry_(false)
      , walkLevel_ (org.walkLevel_)
      , level_(org.level_)
      , geoInFather_ (GeometryImp())
      , geoInFatherImp_(grid_.getRealImplementation(geoInFather_))
      , isLeaf_ (org.isLeaf_)
      , refElem_(grid_.referenceElement())
  {  }

  template<int dim, class GridImp>
  inline void ALU3dGridEntity<0,dim,GridImp> ::
  removeElement ()
  {
    item_  = 0;
  }

  template<int dim, class GridImp>
  inline void ALU3dGridEntity<0,dim,GridImp> ::
  reset (int walkLevel )
  {
    // assert( walkLevel_ >= 0 );

    item_       = 0;
    ghost_      = 0;
    builtgeometry_ = false;
    walkLevel_     = walkLevel;
    level_      = -1;
    isLeaf_     = false;
  }

  // works like assignment
  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp> :: setEntity(const ALU3dGridEntity<0,dim,GridImp> & org)
  {
    item_          = org.item_;
    ghost_         = org.ghost_;
    builtgeometry_ = false;
    level_         = org.level_;
    walkLevel_     = org.walkLevel_;
    isLeaf_        = org.isLeaf_;
  }

  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp>::
  setElement(ALU3DSPACE HElementType & element)
  {
    item_= static_cast<IMPLElementType *> (&element);
    assert( item_ );
    // make sure this method is not called for ghosts
    assert( ! item_->isGhost() );
    ghost_   = 0;
    builtgeometry_=false;
    level_   = (*item_).level();
    isLeaf_  = ((*item_).down() == 0);
  }

  template<int dim, class GridImp>
  inline void
  ALU3dGridEntity<0,dim,GridImp> :: setGhost(ALU3DSPACE HBndSegType & ghost)
  {
    // use element as ghost
    typedef typename ALU3dImplTraits<GridImp::elementType>::IMPLElementType IMPLElementType;
    item_  = static_cast<IMPLElementType *> ( ghost.getGhost().first );

    // method getGhost can return 0, but then is something wrong
    assert(item_);
    assert(item_->isGhost());

    level_   = item_->level();
    // remember pointer to ghost face
    ghost_ = static_cast<PLLBndFaceType *> (&ghost);
    assert( ghost_ );
    builtgeometry_ = false;

    PLLBndFaceType * dwn = static_cast<PLLBndFaceType *> (ghost.down());
    if ( ! dwn ) isLeaf_ = true;
    else
    {
      assert( ghost.level() == level_ );
      if(dwn->ghostLevel() == level_)
        isLeaf_ = true;
      else
        isLeaf_ = false;
    }
    // check wether ghost is leaf or not, ghost leaf means
    // that this is the ghost that we want in the leaf iterator
    // not necessarily is real leaf element
    // see Intersection Iterator, same story
  }

  template<int dim, class GridImp>
  inline int
  ALU3dGridEntity<0,dim,GridImp> :: level() const
  {
    return level_;
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> ::
  equals (const ALU3dGridEntity<0,dim,GridImp> &org ) const
  {
    return (item_ == org.item_);
  }

  template<int dim, class GridImp>
  inline const typename ALU3dGridEntity<0,dim,GridImp>::Geometry &
  ALU3dGridEntity<0,dim,GridImp> :: geometry () const
  {
    assert(item_ != 0);
    if(!builtgeometry_) builtgeometry_ = geoImp_.buildGeom(*item_);
    return geo_;
  }

  template<int dim, class GridImp>
  inline GeometryType
  ALU3dGridEntity<0,dim,GridImp> :: type () const
  {
    return geo_.type();
  }

  template<int dim, class GridImp>
  inline const typename ALU3dGridEntity<0,dim,GridImp>::Geometry &
  ALU3dGridEntity<0,dim,GridImp> :: geometryInFather () const
  {
    assert( item_ );
    const int child = item_->nChild();
    typedef MakeableInterfaceObject<Geometry> GeometryObject;
    typedef typename GeometryObject::ImplementationType GeometryImp;
    // to be improved, when we using not the refine 8 rule
    // see alu3dutility.hh for implementation
    static LocalGeometryStorage<GeometryObject,8> geoms;
    if(!geoms.geomCreated(child))
    {
      typedef typename GridImp::template Codim<0> ::EntityPointer EntityPointer;
      const EntityPointer ep = father();
      geoms.create(grid_,(*ep).geometry(),geometry(),child );
    }
    return geoms[child];
  }

  template<int dim, class GridImp>
  inline int ALU3dGridEntity<0,dim,GridImp> :: getIndex() const
  {
    assert( item_ );
    return (*item_).getIndex();
  }

  //********* begin method subIndex ********************
  // partial specialisation of subIndex
  template <class IMPLElemType, ALU3dGridElementType type, int codim>
  struct IndexWrapper {};

  // specialisation for vertices
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 3>
  {
    typedef ElementTopologyMapping<type> Topo;

    static inline int subIndex(const IMPLElemType &elem, int i)
    {
      return elem.myvertex( Topo::dune2aluVertex(i) )->getIndex(); // element topo
    }
  };

  // specialisation for faces
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type , 1>
  {
    static inline int subIndex(const IMPLElemType &elem, int i)
    {
      // is specialised for each element type and uses
      // the dune2aluFace mapping
      return (getFace(elem,i))->getIndex();
    }
  };

  // specialisation for edges
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 2>
  {
    typedef ElementTopologyMapping<type> Topo;

    // return subIndex of given edge
    static inline int subIndex(const IMPLElemType &elem, int i)
    {
      // get hedge1 corresponding to dune reference element and return number
      return elem.myhedge1( Topo::dune2aluEdge(i) )->getIndex();
    }
  };

  // specialisation for elements
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 0>
  {
    static inline int subIndex(const IMPLElemType &elem, int i) {
      // just return the elements index
      return elem.getIndex();
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  inline int ALU3dGridEntity<0,dim,GridImp> :: getSubIndex (int i) const
  {
    assert(item_ != 0);
    typedef typename  ALU3dImplTraits<GridImp::elementType>::IMPLElementType IMPLElType;
    return IndexWrapper<IMPLElType,GridImp::elementType,cc>::subIndex ( *item_, i);
  }

  //******** end method count *************
  template<int dim, class GridImp>
  template<int cc>
  inline int ALU3dGridEntity<0,dim,GridImp> :: count () const
  {
    return refElem_.size(cc);
  }

  //******** begin method entity ******************
  template <class GridImp, int dim, int cd> struct SubEntities {};

  // specialisation for elements
  template <class GridImp, int dim>
  struct SubEntities<GridImp, dim, 0>
  {
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp>::template Codim<0>::EntityPointer
    entity (const GridImp & grid,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<0, GridImp>(grid , en );
    }
  };

  // specialisation for faces
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,1>
  {
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<1>:: EntityPointer
    entity (const GridImp& grid,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType>::IMPLElementType & item,
            int duneFace)
    {
      int aluFace = Topo::dune2aluFace(duneFace);
      return
        ALU3dGridEntityPointer<1,GridImp>
          (grid,
          level,
          *(getFace(item, duneFace)),    // getFace already constains dune2aluFace
          item.twist(aluFace),
          duneFace    // we need the duneFace number here for the buildGeom method
          );
    }
  };

  // specialisation for edges
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,2>
  {
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;
    typedef typename GridImp::ctype coordType;

    typedef typename GridImp :: ReferenceElementType ReferenceElementType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<2>:: EntityPointer
    entity (const GridImp & grid,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType>::IMPLElementType & item,
            int i)
    {
      // get reference element
      const ReferenceElementType & refElem = grid.referenceElement();

      // get first local vertex number of edge i
      int localNum = refElem.subEntity(i,2,0,dim);

      // get number of first vertex on edge
      int v = en.template getSubIndex<dim> (localNum);

      // get the hedge object
      const typename ALU3dImplTraits<GridImp::elementType>::GEOEdgeType &
      edge = *(item.myhedge1(Topo::dune2aluEdge(i)));

      int vx = edge.myvertex(0)->getIndex();

      // check whether vertex numbers are equal, otherwise twist is 1
      int twst = (v != vx) ? 1 : 0;
      return ALU3dGridEntityPointer<2,GridImp> (grid, level, edge, twst );
    }
  };

  // specialisation for vertices
  template <class GridImp, int dim>
  struct SubEntities<GridImp,dim,3>
  {
    typedef ElementTopologyMapping<GridImp::elementType> Topo;
    typedef ALU3dGridEntity<0,dim,GridImp> EntityType;

    static typename ALU3dGridEntity<0,dim,GridImp> :: template Codim<3>:: EntityPointer
    entity (const GridImp & grid,
            const int level,
            const EntityType & en,
            const typename ALU3dImplTraits<GridImp::elementType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<3,GridImp>
               (grid, level, *item.myvertex(Topo::dune2aluVertex(i))); // element topo
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  inline typename ALU3dGridEntity<0,dim,GridImp>::template Codim<cc>:: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: entity (int i) const
  {
    return SubEntities<GridImp,dim,cc>::entity(grid_,level(),*this,*item_,i);
  }

  //**** end method entity *********

  template<int dim, class GridImp>
  inline PartitionType ALU3dGridEntity<0,dim,GridImp> ::
  partitionType () const
  {
    assert( item_ );
    // make sure we really got a ghost
    assert( (isGhost()) ? item_->isGhost() : true );
    return (isGhost() ?  GhostEntity : InteriorEntity);
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: isLeaf() const
  {
    return isLeaf_;
  }
  template<int dim, class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp>
  ALU3dGridEntity<0,dim,GridImp> :: hbegin (int maxlevel) const
  {
    assert(item_ != 0);
    // if isGhost is true the end iterator will be returned
    return ALU3dGridHierarchicIterator<GridImp>(grid_,*item_,maxlevel, isGhost() );
  }

  template<int dim, class GridImp>
  inline ALU3dGridHierarchicIterator<GridImp> ALU3dGridEntity<0,dim,GridImp> :: hend (int maxlevel) const
  {
    assert(item_ != 0);
    return ALU3dGridHierarchicIterator<GridImp> (grid_,*item_,maxlevel,true);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLeafIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ileafbegin () const
  {
    assert(item_ != 0);
    // NOTE: normaly here false should be given, which means that we create a non
    // end iterator, but isGhost() is normaly false. If isGhost() is true,
    // an end iterator is created,
    // because on ghosts we dont run itersection iterators
    return ALU3dGridIntersectionIteratorType (grid_,*this,walkLevel_, isGhost() );
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLeafIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ileafend () const
  {
    assert(item_ != 0);
    return ALU3dGridLeafIntersectionIteratorType (grid_, *this ,walkLevel_,true);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLevelIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ilevelbegin () const
  {
    assert(item_ != 0);
    // NOTE: normaly here false should be given, which means that we create a non
    // end iterator, but isGhost() is normaly false. If isGhost() is true,
    // an end iterator is created,
    // because on ghosts we dont run itersection iterators
    return ALU3dGridLevelIntersectionIteratorType (grid_,*this,walkLevel_, isGhost() );
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: ALU3dGridLevelIntersectionIteratorType
  ALU3dGridEntity<0,dim,GridImp> :: ilevelend () const
  {
    assert(item_ != 0);
    return ALU3dGridLevelIntersectionIteratorType (grid_, *this ,walkLevel_,true);
  }

  template<int dim, class GridImp>
  inline typename ALU3dGridEntity<0,dim,GridImp> :: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: father() const
  {
    if(! item_->up() )
    {
      std::cerr << "ALU3dGridEntity<0," << dim << "," << dimworld << "> :: father() : no father of entity globalid = " << getIndex() << "\n";
      return ALU3dGridEntityPointer<0,GridImp> (grid_, static_cast<ALU3DSPACE HElementType &> (*item_));
    }
    return ALU3dGridEntityPointer<0,GridImp> (grid_, static_cast<ALU3DSPACE HElementType &> (*(item_->up())));
  }

  // Adaptation methods
  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: mark (int ref) const
  {
    assert(item_ != 0);

    // if this assertion is thrown then you try to mark a non leaf entity
    // which is leads to unpredictable results
    if( !isLeaf() ) return false;

    // mark for coarsening
    if(ref < 0)
    {
      // don't mark macro elements for coarsening ;)
      if(level() <= 0) return false;

      (*item_).request(coarse_element_t);
      return true;
    }

    // mark for refinement
    if(ref > 0)
    {
      (*item_).request(refine_element_t);
      return true;
    }

    // mark for none
    (*item_).request( nosplit_element_t );
    return true;
  }

  // return mark of entity
  template<int dim, class GridImp>
  inline int ALU3dGridEntity<0,dim,GridImp> :: getMark () const
  {
    assert(item_ != 0);

    const MarkRuleType rule = (*item_).requestrule();

    if(rule == refine_element_t) return 1;
    if(rule == coarse_element_t) return -1;
    assert( rule == nosplit_element_t );
    return 0;
  }

  // Adaptation methods
  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: isNew () const
  {
    assert( item_ );
    return item_->hasBeenRefined();
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: mightVanish () const
  {
    assert( item_ );
    return ((*item_).requestrule() == coarse_element_t);
  }

  template<int dim, class GridImp>
  inline bool ALU3dGridEntity<0,dim,GridImp> :: hasBoundaryIntersections () const
  {
    // on ghost elements return false
    if( isGhost() ) return false;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    typedef typename ALU3dImplTraits<GridImp::elementType>::HasFaceType HasFaceType;
    typedef typename ALU3dImplTraits<GridImp::elementType>::GEOFaceType GEOFaceType;

    assert( item_ );
    for(int i=0; i<numFaces; ++i)
    {
      const GEOFaceType & face = *getFace(*item_,i);

#if ALU3DGRID_PARALLEL
      // don't count internal boundaries as boundary
      if( face.isBorder() ) continue ;
#endif

      // check both
      const HasFaceType * outerElement = face.nb.front().first;
      // if we got our own element, get other side
      if( item_ == outerElement )
      {
        outerElement = face.nb.rear().first;
      }

      assert( outerElement );
      if( outerElement->isboundary() ) return true;
    }
    return false;
  }

  //*******************************************************************
  //
  //  --EntityPointer
  //  --EnPointer
  //
  //*******************************************************************

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid,
                             const int level,
                             const MyHElementType &item)
    : grid_(grid)
      , item_(const_cast<MyHElementType *> (&item))
      , entity_(0)
      , locked_(false)
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid,
                             const HBndSegType & ghostFace )
    : grid_(grid)
      , item_(0)
      , entity_ ( grid_.template getNewEntity<codim> ( ghostFace.level() ))
      , locked_( true ) // entity should not be released
  {
    // sets entity and item pointer
    updateGhostPointer( const_cast<HBndSegType &> (ghostFace) );
  }

  // constructor Level,Leaf and HierarchicIterator
  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid, int level )
    : grid_(grid)
      , item_(0)
      , entity_ ( grid_.template getNewEntity<codim> ( level ) )
      , locked_(true) // entity should not be released
  {
    // this needs to be called
    // have to investigate why
    entityImp().reset(level);
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const ALU3dGridEntityPointerType & org)
    : grid_(org.grid_)
      , item_(org.item_)
      , entity_(0)
      , locked_(org.locked_)
  {
    // if entity exists then copy entity
    getEntity( org );
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp> ::
  getEntity(const ALU3dGridEntityPointerType & org)
  {
    if( org.entity_ )
    {
      assert( entity_ == 0 );
      entity_ = grid_.template getNewEntity<codim> ();
      this->entityImp().setEntity( org.entityImp() );
    }
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> &
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    assert( &grid_ == &org.grid_ );
    // if entity exists, just free and reset pointers
    if(entity_) this->done();
    // set item
    item_ = org.item_;
    return *this;
  }

  template<int codim, class GridImp >
  inline void
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    assert( &grid_ == &org.grid_ );

    // set item
    item_ = org.item_;

    if(item_)
    {
      // if no entity check org entity
      // if no org entity then nothing is done
      if( !entity_ )
      {
        getEntity(org);
      }
      else
      {
#if ALU3DGRID_PARALLEL
        // in case of ghost element use different set method
        if( item_->isGhost() )
        {
          // on ghosts entity pointers entity always exists
          assert( org.entity_ );
          this->entityImp().setEntity( org.entityImp() );
        }
        else
#endif
        {
          // otherwise item is set
          this->entityImp().setElement( *item_ );
        }
      }
    }
    else
    {
      this->done();
    }
    return ;
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointerBase<codim,GridImp> ::
  ~ALU3dGridEntityPointerBase()
  {
    this->done();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::done ()
  {
    item_ = 0;
    // free entity
    freeEntity();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::freeEntity ()
  {
    // sets entity pointer in the status of an empty entity
    if( entity_ )
    {
      this->entityImp().removeElement();
      grid_.template freeEntity<codim> ( (EntityObject *) entity_ );
      entity_ = 0;
    }
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::compactify()
  {
    // sets entity pointer in the status of an empty entity
    if( ! locked_ )
    {
      freeEntity ();
    }
  }

  template<int codim, class GridImp >
  inline bool ALU3dGridEntityPointerBase<codim,GridImp>::
  equals (const ALU3dGridEntityPointerBase<codim,GridImp>& i) const
  {
    // check equality of underlying items
    return (item_ == i.item_);
  }

  template<int codim, class GridImp >
  inline typename ALU3dGridEntityPointerBase<codim,GridImp>::Entity &
  ALU3dGridEntityPointerBase<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( item_ );
    if(!entity_)
    {
      entity_ = grid_.template getNewEntity<codim> ();
      this->entityImp().setElement( *item_ );
    }
    assert( item_ == & this->entityImp().getItem() );
    return (*entity_);
  }

  template<int codim, class GridImp >
  inline int ALU3dGridEntityPointerBase<codim,GridImp>::level () const
  {
    assert( item_ );
    return item_->level();
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateGhostPointer( ALU3DSPACE HBndSegType & ghostFace )
  {
    assert( entity_ );
    this->entityImp().setGhost( ghostFace );
    // inside the method setGhost the method getGhost of the ghostFace is
    // called and set as item
    const MyHElementType * item = & (this->entityImp().getItem());
    item_ = const_cast<MyHElementType *> (item);
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateEntityPointer( MyHElementType * item , int )
  {
    item_ = item;
    if( item_ && entity_ )
    {
      this->entityImp().setElement( *item_ );
    }
  }

  ///////////////////////////////////////////////////////////////////
  //
  //  specialisation for higher codims
  //
  ///////////////////////////////////////////////////////////////////

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const GridImp & grid,
                         const int level,
                         const MyHElementType &item,
                         const int twist,
                         const int duneFace )
    : ALU3dGridEntityPointerBase<codim,GridImp> (grid,level,item)
      , level_(level)
      , twist_ (twist)
      , face_(duneFace)
  {
    assert( (codim == 1) ? (face_ >= 0) : 1 );
  }

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const ALU3dGridEntityPointerType & org)
    : ALU3dGridEntityPointerBase<codim,GridImp>(org)
      , level_(org.level_)
      , twist_(org.twist_)
      , face_(org.face_)
  {}

  template<int codim, class GridImp >
  inline ALU3dGridEntityPointer<codim,GridImp> &
  ALU3dGridEntityPointer<codim,GridImp>::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    // first copy level, twist and face, because this might be used in
    // clone
    level_ = org.level_;
    twist_ = org.twist_;
    face_  = org.face_;

    // copy item pointer
    clone(org);
    return *this;
  }

  template<int codim, class GridImp >
  inline void
  ALU3dGridEntityPointer<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    assert( &this->grid_ == &org.grid_ );
    // set item
    this->item_ = org.item_;

    // if entity exists, just remove item pointer
    if(this->item_)
    {
      if( ! this->entity_ )
        this->getEntity(org);
      else
        this->entityImp().setElement( *this->item_ , this->level(), twist_ , face_ );
    }
    else
      this->done();
    return ;
  }

  template<int codim, class GridImp >
  inline typename ALU3dGridEntityPointer<codim,GridImp>::Entity &
  ALU3dGridEntityPointer<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( this->item_ );
    if(!this->entity_)
    {
      this->entity_ = this->grid_.template getNewEntity<codim> ();
      this->entityImp().setElement( *this->item_ , this->level(), twist_ , face_ );
    }
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }

  template<int codim, class GridImp >
  inline int ALU3dGridEntityPointer<codim,GridImp>::level () const
  {
    return level_;
  }

  template<int codim, class GridImp >
  inline void ALU3dGridEntityPointer<codim,GridImp>::
  updateEntityPointer( MyHElementType * item, int level)
  {
    this->item_ = item;
    level_ = level;
    if( this->item_ && this->entity_ )
    {
      this->entityImp().setElement( *this->item_ , level_ );
    }
  }


} // end namespace Dune
