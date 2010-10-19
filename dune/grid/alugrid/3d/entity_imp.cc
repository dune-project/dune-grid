// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_ENTITY_CC
#define DUNE_ALUGRID_ENTITY_CC

#if COMPILE_ALUGRID_INLINE == 0
#include <config.h>
#endif

#include "alu3dinclude.hh"
#include <dune/grid/alugrid/3d/alugrid.hh>
#include "entity.hh"

#if COMPILE_ALUGRID_INLINE == 0
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif
#include <dune/grid/alugrid/common/geostorage.hh>

#if COMPILE_ALUGRID_INLINE
#define alu_inline inline
#else
#define alu_inline
#endif

namespace Dune {

  // --Entity
  template <int cd, int dim, class GridImp>
  ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const GridImp  &grid, int level) :
    geo_( GeometryImp() ),
    grid_(grid),
    item_(0),
    level_(0),
    gIndex_(-1),
    twist_(0),
    face_(-1),
    builtgeometry_(false),
    partitionType_(InteriorEntity)
  {}

  // --Entity
  template <int cd, int dim, class GridImp>
  alu_inline ALU3dGridEntity<cd,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<cd,dim,GridImp> & org) :
    geo_( GeometryImp() ),
    grid_(org.grid_),
    item_(org.item_),
    level_(org.level_),
    gIndex_(org.gIndex_),
    twist_(org.twist_),
    face_(org.face_),
    builtgeometry_(false),
    partitionType_(org.partitionType_)
  {}

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
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
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const HItemType & item)
  {
    setElement(item,GetLevel<GridImp,cd>::getLevel(grid_,item));
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const ALU3dGridEntityKeyType& key)
  {
    setElement(*key.item(), key.level(), key.twist(), key.face());
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setElement(const HItemType & item, const int level, int twist , int face )
  {
    // cast interface to implementation
    item_   = static_cast<const ItemType *> (&item);
    gIndex_ = (*item_).getIndex();
    twist_  = twist;
    level_  = level;
    face_   = face;
    builtgeometry_=false;
    partitionType_ = this->convertBndId( *item_ );
  }

  template<int cd, int dim, class GridImp>
  alu_inline void ALU3dGridEntity<cd,dim,GridImp> ::
  setGhost(const HBndSegType &ghost)
  {
    // this method only exists, that we don't have to specialize the
    // Iterators for each codim, this method should not be called otherwise
    // error
    DUNE_THROW(GridError,"This method should not be called!");
  }

  template<int cd, int dim, class GridImp>
  alu_inline PartitionType ALU3dGridEntity<cd,dim,GridImp> ::
  convertBndId(const HItemType & item) const
  {
    if(item.isGhost())
    {
      return GhostEntity;
    }
    else if(item.isBorder())
    {
      return BorderEntity;
    }
    else
    {
      assert( item.isInterior() );
      return InteriorEntity;
    }
  }

  template<int cd, int dim, class GridImp>
  alu_inline const typename ALU3dGridEntity<cd,dim,GridImp>::Geometry &
  ALU3dGridEntity<cd,dim,GridImp>:: geometry() const
  {
    //assert( (cd == 1) ? (face_ >= 0) : 1 );
    if(!builtgeometry_) builtgeometry_ = geoImp().buildGeom(*item_, twist_, face_ );
    return geo_;
  }

  /////////////////////////////////////////////////
  //
  //  --Entity0
  //  --Codim0Entity
  //
  /////////////////////////////////////////////////

  template<int dim, class GridImp>
  alu_inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const GridImp  &grid, int wLevel)
    : geo_(GeometryImp())
      , grid_( grid )
      , item_( 0 )
      , ghost_( 0 )
      , level_(-1)
      , isLeaf_ (false)
      , builtgeometry_(false)
  {  }

  template<int dim, class GridImp>
  alu_inline ALU3dGridEntity<0,dim,GridImp> ::
  ALU3dGridEntity(const ALU3dGridEntity<0,dim,GridImp> & org)
    : geo_(GeometryImp())
      , grid_(org.grid_)
      , item_(org.item_)
      , ghost_( org.ghost_ )
      , level_(org.level_)
      , isLeaf_ (org.isLeaf_)
      , builtgeometry_(false)
  {  }

  template<int dim, class GridImp>
  alu_inline const typename ALU3dGridEntity<0,dim,GridImp>::Geometry &
  ALU3dGridEntity<0,dim,GridImp> :: geometry () const
  {
    assert(item_ != 0);
    if( ! builtgeometry_ ) builtgeometry_ = geoImp().buildGeom(*item_);
    return geo_;
  }

  template<int dim, class GridImp>
  alu_inline const typename ALU3dGridEntity<0,dim,GridImp>::Geometry &
  ALU3dGridEntity<0,dim,GridImp> :: geometryInFather () const
  {
    assert( item_ );
    const int child = item_->nChild();
    typedef MakeableInterfaceObject<Geometry> GeometryObject;
    typedef typename GeometryObject::ImplementationType GeometryImp;
    // to be improved, when we using not the refine 8 rule
    // see dune/grid/alugrid/common/geostrage.hh for implementation
    static ALULocalGeometryStorage<GridImp, GeometryObject, 8> geoms( type(), true);
    assert( geoms.geomCreated(child) );
    return geoms[ child ];
  }

  //********* begin method subIndex ********************
  // partial specialisation of subIndex
  template <class IMPLElemType, ALU3dGridElementType type, int codim>
  struct IndexWrapper {};

  // specialisation for vertices
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 3>
  {
    typedef ElementTopologyMapping<type> ElemTopo;

    static int subIndex(const IMPLElemType &elem, int i)
    {
      return elem.myvertex( ElemTopo::dune2aluVertex(i) )->getIndex(); // element topo
    }
  };

  // specialisation for faces
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type , 1>
  {
    static int subIndex(const IMPLElemType &elem, int i)
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
    typedef ElementTopologyMapping<type> ElemTopo;

    // return subIndex of given edge
    static int subIndex(const IMPLElemType &elem, int i)
    {
      // get hedge1 corresponding to dune reference element and return number
      return elem.myhedge1( ElemTopo::dune2aluEdge(i) )->getIndex();
    }
  };

  // specialisation for elements
  template <class IMPLElemType, ALU3dGridElementType type>
  struct IndexWrapper<IMPLElemType, type, 0>
  {
    static int subIndex(const IMPLElemType &elem, int i) {
      // just return the elements index
      return elem.getIndex();
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: getSubIndex (int i) const
  {
    assert(item_ != 0);
    typedef typename  ImplTraits::IMPLElementType IMPLElType;
    return IndexWrapper<IMPLElType,GridImp::elementType,cc>::subIndex ( *item_, i);
  }

  template<int dim, class GridImp>
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: subIndex (int i, unsigned int codim ) const
  {
    typedef ElementTopologyMapping<GridImp::elementType> ElemTopo;

    assert(item_ != 0);
    switch (codim)
    {
    case 0 :
      return this->getIndex();
    case 1 :
      return (ALU3dGridFaceGetter< Comm >::getFace( *item_, i ))->getIndex();
    case 2 :
      return item_->myhedge1( ElemTopo::dune2aluEdge( i ) )->getIndex();
    case 3 :
      return item_->myvertex( ElemTopo::dune2aluVertex( i ) )->getIndex();
    default :
      assert(false);
      abort();
    }
    return -1;
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
            const EntityType & entity,
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<0, GridImp>( entity );
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
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int duneFace)
    {
      int aluFace = Topo::dune2aluFace(duneFace);
      return
        ALU3dGridEntityPointer<1,GridImp>
          (grid,
          level,
          *(ALU3dGridFaceGetter< typename GridImp::MPICommunicatorType >::getFace(item, duneFace)),    // getFace already constains dune2aluFace
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
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      // get reference element
      const ReferenceElementType & refElem = grid.referenceElement();

      // get first local vertex number of edge i
      int localNum = refElem.subEntity(i,2,0,dim);

      // get number of first vertex on edge
      int v = en.template getSubIndex<dim> (localNum);

      // get the hedge object
      const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::GEOEdgeType &
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
            const typename ALU3dImplTraits<GridImp::elementType, typename GridImp::MPICommunicatorType>::IMPLElementType & item,
            int i)
    {
      return ALU3dGridEntityPointer<3,GridImp>
               (grid, level, *item.myvertex(Topo::dune2aluVertex(i))); // element topo
    }
  };

  template<int dim, class GridImp>
  template<int cc>
  typename ALU3dGridEntity<0,dim,GridImp>::template Codim<cc>:: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: subEntity (int i) const
  {
    return SubEntities<GridImp,dim,cc>::entity(grid_,level(),*this,*item_,i);
  }
  //**** end method entity *********

  template<int dim, class GridImp>
  typename ALU3dGridEntity<0,dim,GridImp> :: EntityPointer
  ALU3dGridEntity<0,dim,GridImp> :: father() const
  {
    HElementType* up = item_->up();
    if( ! up )
    {
      std::cerr << "ALU3dGridEntity<0," << dim << "," << dimworld << "> :: father() : no father of entity globalid = " << getIndex() << "\n";
      return ALU3dGridEntityPointer<0,GridImp> (grid_, static_cast<HElementType &> (*item_));
    }

    if( isGhost () )
    {
      return ALU3dGridEntityPointer<0,GridImp> (grid_, static_cast<const HBndSegType &> (*(getGhost().up())));
    }

    return ALU3dGridEntityPointer<0,GridImp> (grid_, static_cast<HElementType &> ( *up ));
  }

  // Adaptation methods
  template<int dim, class GridImp>
  bool ALU3dGridEntity<0,dim,GridImp> :: mark (int ref) const
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
  alu_inline int ALU3dGridEntity<0,dim,GridImp> :: getMark () const
  {
    assert(item_ != 0);

    const MarkRuleType rule = (*item_).requestrule();

    if(rule == refine_element_t) return 1;
    if(rule == coarse_element_t) return -1;
    assert( rule == nosplit_element_t );
    return 0;
  }


  template<int dim, class GridImp>
  bool ALU3dGridEntity<0,dim,GridImp> :: hasBoundaryIntersections () const
  {
    // on ghost elements return false
    if( isGhost() ) return false;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    typedef typename ImplTraits::HasFaceType HasFaceType;
    typedef typename ImplTraits::GEOFaceType GEOFaceType;

    assert( item_ );
    for(int i=0; i<numFaces; ++i)
    {
      const GEOFaceType &face = *ALU3dGridFaceGetter< Comm >::getFace( *item_, i );

      // don't count internal boundaries as boundary
      if( face.isBorder() ) continue ;

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
#if 0
  template<int codim, class GridImp >
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid,
                             const HElementType &item)
    : grid_(grid)
      , item_(const_cast<HElementType *> (&item))
      , entity_( 0 )
      , locked_ ( false ) // entity can be released
  {}

  template<int codim, class GridImp >
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid,
                             const HBndSegType & ghostFace )
    : grid_(grid)
      , item_(0)
      , entity_ ( grid_.template getNewEntity<codim> ( ghostFace.level() ))
      , locked_( true ) // entity should not be released, otherwise is ghost info lost
  {
    // sets entity and item pointer
    updateGhostPointer( const_cast<HBndSegType &> (ghostFace) );
  }

  // constructor Level,Leaf and HierarchicIterator
  template<int codim, class GridImp >
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const GridImp & grid, int level )
    : grid_(grid)
      , item_(0)
      , entity_ ( grid_.template getNewEntity<codim> ( level ) )
      , locked_ ( false ) // entity can be released
  {
    // this needs to be called
    // have to investigate why
    entityImp().reset(level);
  }

  template<int codim, class GridImp >
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  ALU3dGridEntityPointerBase(const ALU3dGridEntityPointerType & org)
    : grid_(org.grid_)
      , item_(org.item_)
      , entity_( 0 )
      , locked_( org.locked_ )
  {
    // if entity exists then copy entity
    getEntity( org );
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp> ::
  getEntity(const ALU3dGridEntityPointerType & org)
  {
    // if entity existed for original pointer then copy
    if( org.entity_ )
    {
      assert( entity_ == 0 );
      entity_ = grid_.template getNewEntity<codim> ();
      // set entity right away
      this->entityImp().setEntity( org.entityImp() );
    }
  }

  template<int codim, class GridImp >
  ALU3dGridEntityPointerBase<codim,GridImp> &
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    clone( org );
    return *this;
  }

  template<int codim, class GridImp >
  alu_inline void
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    assert( &grid_ == &org.grid_ );

    // set item
    item_   = org.item_;
    // copy locked info
    locked_ = org.locked_;

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
        // in case of ghost element use different set method
        if( item_->isGhost() )
        {
          // on ghosts entity pointers entity always exists
          assert( org.entity_ );
          this->entityImp().setEntity( org.entityImp() );
          locked_ = true ;
        }
        else
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
  ALU3dGridEntityPointerBase<codim,GridImp> ::
  ~ALU3dGridEntityPointerBase()
  {
    this->done();
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp>::done ()
  {
    item_   = 0;
    locked_ = false;
    // free entity
    freeEntity();
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp>::freeEntity ()
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
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp>::compactify()
  {
    // sets entity pointer in the status of an empty entity
    if( ! locked_ )
    {
      freeEntity ();
    }
  }

  template<int codim, class GridImp >
  bool ALU3dGridEntityPointerBase<codim,GridImp>::
  equals (const ALU3dGridEntityPointerBase<codim,GridImp>& i) const
  {
    // check equality of underlying items
    return (item_ == i.item_);
  }

  template<int codim, class GridImp >
  typename ALU3dGridEntityPointerBase<codim,GridImp>::Entity &
  ALU3dGridEntityPointerBase<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( item_ );
    assert( (item_->isGhost()) ? locked_ : true );
    assert( (locked_) ? (entity_ != 0) : true);
    if( ! entity_ )
    {
      entity_ = grid_.template getNewEntity<codim> ();
      this->entityImp().setElement( *item_ );
    }
    assert( item_ == & this->entityImp().getItem() );
    return (*entity_);
  }

  template<int codim, class GridImp >
  alu_inline int ALU3dGridEntityPointerBase<codim,GridImp>::level () const
  {
    assert( item_ );
    return item_->level();
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateGhostPointer( HBndSegType & ghostFace )
  {
    assert( entity_ );
    this->entityImp().setGhost( ghostFace );
    // inside the method setGhost the method getGhost of the ghostFace is
    // called and set as item
    const HElementType * item = & (this->entityImp().getItem());
    item_ = const_cast<HElementType *> (item);
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointerBase<codim,GridImp>::
  updateEntityPointer( HElementType * item , int )
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
  alu_inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const GridImp & grid,
                         const int level,
                         const HElementType &item,
                         const int twist,
                         const int duneFace )
    : ALU3dGridEntityPointerBase<codim,GridImp> (grid,item)
      , level_(level)
      , twist_ (twist)
      , face_(duneFace)
  {
    assert( (codim == 1) ? (face_ >= 0) : 1 );
  }

  template<int codim, class GridImp >
  alu_inline ALU3dGridEntityPointer<codim,GridImp> ::
  ALU3dGridEntityPointer(const ALU3dGridEntityPointerType & org)
    : ALU3dGridEntityPointerBase<codim,GridImp>(org)
      , level_(org.level_)
      , twist_(org.twist_)
      , face_(org.face_)
  {}

  template<int codim, class GridImp >
  alu_inline ALU3dGridEntityPointer<codim,GridImp> &
  ALU3dGridEntityPointer<codim,GridImp>::
  operator = (const ALU3dGridEntityPointerType & org)
  {
    // clone pointer
    clone(org);
    return *this;
  }

  template<int codim, class GridImp >
  alu_inline void
  ALU3dGridEntityPointer<codim,GridImp> ::
  clone (const ALU3dGridEntityPointerType & org)
  {
    // first copy level, twist and face, because this might be used in
    // clone
    level_ = org.level_;
    twist_ = org.twist_;
    face_  = org.face_;

    assert( &this->grid_ == &org.grid_ );
    // set item
    this->item_ = org.item_;
    // copy lock status
    this->locked_ = org.locked_;

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
  alu_inline typename ALU3dGridEntityPointer<codim,GridImp>::Entity &
  ALU3dGridEntityPointer<codim,GridImp>::dereference () const
  {
    // don't dereference empty entity pointer
    assert( this->item_ );
    if( ! this->entity_ )
    {
      this->entity_ = this->grid_.template getNewEntity<codim> ();
      this->entityImp().setElement( *this->item_ , this->level(), twist_ , face_ );
    }
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }

  template<int codim, class GridImp >
  alu_inline int ALU3dGridEntityPointer<codim,GridImp>::level () const
  {
    return level_;
  }

  template<int codim, class GridImp >
  alu_inline void ALU3dGridEntityPointer<codim,GridImp>::
  updateEntityPointer( HElementType * item, int level)
  {
    this->item_ = item;
    level_ = level;
    if( this->item_ && this->entity_ )
    {
      this->entityImp().setElement( *this->item_ , level_ );
    }
  }
#endif

#if COMPILE_ALUGRID_LIB
  // Instantiation
  template class ALU3dGrid< hexa, No_Comm >;
  template class ALU3dGrid< tetra, No_Comm >;

  // Instantiation
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<1, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<1, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<2, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<2, 3, const ALU3dGrid< hexa, No_Comm > >;

  template class ALU3dGridEntity<3, 3, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridEntity<3, 3, const ALU3dGrid< hexa, No_Comm > >;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 0 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 0 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 1 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 1 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 2 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 2 >( int ) const;

  template ALU3dGrid< tetra, No_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, No_Comm > > :: subEntity< 3 >( int ) const;
  template ALU3dGrid< hexa, No_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, No_Comm > > :: subEntity< 3 >( int ) const;

#if ALU3DGRID_PARALLEL
  // Instantiation
  template class ALU3dGrid< hexa, MPI_Comm >;
  template class ALU3dGrid< tetra, MPI_Comm >;

  // Instantiation with MPI
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<1, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<1, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<2, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<2, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template class ALU3dGridEntity<3, 3, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridEntity<3, 3, const ALU3dGrid< hexa, MPI_Comm > >;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 0 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 0 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 0 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 1 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 1 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 1 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 2 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 2 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 2 >( int ) const;

  template ALU3dGrid< tetra, MPI_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< tetra, MPI_Comm > > :: subEntity< 3 >( int ) const;
  template ALU3dGrid< hexa, MPI_Comm > :: Traits :: Codim< 3 > :: EntityPointer
  ALU3dGridEntity<0, 3, const ALU3dGrid< hexa, MPI_Comm > > :: subEntity< 3 >( int ) const;
#endif // #if ALU3DGRID_PARALLEL

#endif // #if COMPILE_ALUGRID_LIB
}
#undef alu_inline
#endif
