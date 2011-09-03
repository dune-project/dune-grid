// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDENTITYIMP_CC
#define DUNE_ALU2DGRIDENTITYIMP_CC

#include "geometry.hh"
#include "grid.hh"
#include <dune/common/exceptions.hh>

#include <dune/grid/alugrid/common/geostorage.hh>

namespace Dune {


  //**********************************************************************
  //
  // --ALU2dGridEntity
  // --Entity
  //
  //**********************************************************************


  //! level of this element
  template<int cd, int dim, class GridImp>
  inline int ALU2dGridEntity<cd, dim, GridImp> ::
  level () const {
    if (level_ == -1) {
      if (cd != 2)
        level_ = item_->level();
      else
        level_ = item_->level() + 1;
    }
    assert(level_ != -1);
    return level_;
  }

  // forward declararion of struct ElementWrapper
  template <int codim, int dim, class GridImp>
  struct ElementWrapper;


  template<int cd, int dim, class GridImp>
  inline bool ALU2dGridEntity<cd, dim, GridImp> :: equals(const ALU2dGridEntity<cd, dim, GridImp> &org) const {
    return ElementWrapper<cd,dim, GridImp>::isTheSame (*item_, face_, *org.item_, org.face_);
  }

  //! Constructor
  template<int cd, int dim, class GridImp>
  inline ALU2dGridEntity<cd, dim, GridImp> ::
  ALU2dGridEntity(const FactoryType& factory, int level)
    : factory_( factory ),
      item_(0),
      geoObj_(GeometryImp()),
      level_(level),
      face_(-1)
  {}


  template<int cd, int dim, class GridImp>
  inline void ALU2dGridEntity<cd,dim,GridImp>:: setElement(const ElementType &element, int face, int level) const {
    item_= const_cast<ElementType *> (&element);
    level_ = level;
    face_ = face;

    geoImpl().unsetUp2Date();
  }


  template<int cd, int dim, class GridImp>
  inline void ALU2dGridEntity<cd,dim,GridImp>:: setElement( const EntitySeed& seed ) const
  {
    setElement( *(seed.item()), seed.face(), seed.level() );
  }

  //! set item pointer to NULL
  template<int cd, int dim, class GridImp>
  inline void ALU2dGridEntity<cd,dim,GridImp> :: removeElement() {
    item_ = 0;
  }

  //! Copy Constructor
  template<int cd, int dim, class GridImp>
  inline ALU2dGridEntity<cd, dim, GridImp> ::
  ALU2dGridEntity(const ALU2dGridEntity<cd,dim,GridImp> & org)
    : factory_( org.factory_ ),
      item_(org.item_),
      geoObj_(GeometryImp()),
      level_(org.level_),
      face_(org.face_)
  {}

  //! geometry of this entity
  template<int cd, int dim, class GridImp>
  inline const typename ALU2dGridEntity<cd, dim, GridImp> :: Geometry &
  ALU2dGridEntity<cd, dim, GridImp> :: geometry () const
  {
    if( !geoImpl().up2Date() )
      geoImpl().buildGeom(*item_,face_);

    assert( geoImpl().up2Date() );
    return geoObj_;
  }

  //! geometry type of geometry of this entity
  template<int cd, int dim, class GridImp>
  inline GeometryType
  ALU2dGridEntity<cd, dim, GridImp> :: type () const
  {
    return geoObj_.type();
  }

  template<int cd, int dim, class GridImp>
  inline int ALU2dGridEntity<cd,dim,GridImp > :: getIndex() const
  {
    assert(item_ != 0);
    return ElementWrapper<cd, dim, GridImp>::getElemIndex (grid(), *item_, face_);
  }

  /**
     \brief Id of the boundary which is associated with
     the entity, returns 0 for inner entities, arbitrary int otherwise
   */

  template<int cd, int dim, class GridImp>
  inline int ALU2dGridEntity<cd,dim,GridImp> :: boundaryId() const {
    int isBoundary=0, i=0;
    while(!isBoundary && i<dim) {
      if (item_->nbbnd(i)!=0)
        isBoundary = item_->nbbnd(i)->type();
      ++i;
    }
    return isBoundary;
  }

  //**********************************************************************
  //
  //  --ALU2dGridEntity
  //  --0Entity
  //
  //**********************************************************************

  //! Constructor creating empty Entity
  template<int dim, class GridImp>
  inline ALU2dGridEntity<0,dim,GridImp> ::
  ALU2dGridEntity(const FactoryType& factory, int level)
    : factory_( factory )
      , item_(0)
      , geoObj_(GeometryImp())
      , isLeaf_ (false)
  {}

  //! Copy Constructor
  template<int dim, class GridImp>
  inline ALU2dGridEntity<0, dim, GridImp> ::
  ALU2dGridEntity(const ALU2dGridEntity<0,dim,GridImp> & org)
    : factory_( org.factory_ ),
      item_(org.item_),
      geoObj_(GeometryImp()),
      isLeaf_(org.isLeaf_)
  {}

  //! level of this element
  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: level () const {
    assert( item_ );
    return (*item_).level();
  }

  //! geometry of this entity
  template<int dim, class GridImp>
  inline const typename ALU2dGridEntity<0, dim, GridImp> :: Geometry & ALU2dGridEntity<0,dim,GridImp> ::
  geometry () const
  {
    assert(item_ != 0);
    if(! geoImpl().up2Date() )
      geoImpl().buildGeom(*item_);

    assert( geoImpl().up2Date() );
    return geoObj_;
  }

  //! geometry type of geometry of this entity
  template<int dim, class GridImp>
  inline GeometryType
  ALU2dGridEntity<0, dim, GridImp> :: type () const
  {
    return geoObj_.type();
  }

  //! returns true if Entity is leaf (i.e. has no children)
  template<int dim, class GridImp>
  inline bool ALU2dGridEntity<0,dim,GridImp> :: isLeaf () const {
    return isLeaf_;
  }

  //! Inter-level access to father element on coarser grid.
  //! Assumes that meshes are nested.
  template<int dim, class GridImp>
  inline typename ALU2dGridEntity<0, dim, GridImp> :: EntityPointer ALU2dGridEntity<0,dim,GridImp> :: father () const {
    // don't request for father on macro level
    assert(level()>0);
    return EntityPointer(factory_, *(item_->father()));
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0, dim, GridImp> :: nChild() const
  {
    assert( item_ );
    return item_->childNr();
  }

  template<int dim, class GridImp>
  inline const typename ALU2dGridEntity<0, dim, GridImp>:: LocalGeometry & ALU2dGridEntity<0,dim,GridImp> ::
  geometryInFather () const
  {
    assert( level() > 0 );

    typedef MakeableInterfaceObject<LocalGeometry> LocalGeometryObject;

    const GeometryType myType = type();
    // we need to storages in case of cube grid,
    // one for quadrilaterals and one for triangles
    if( GridImp :: elementType != ALU2DSPACE triangle && myType.isCube() )
    {
      assert( grid().nonConform() );
      static ALULocalGeometryStorage<GridImp, LocalGeometryObject,4> geoms( myType, true );
      return geoms[ nChild() ];
    }
    else
    {
      if( grid().nonConform() )
      {
        static ALULocalGeometryStorage<GridImp, LocalGeometryObject,4> geoms( myType, true );
        return geoms[ nChild() ];
      }
      else
      {
        static ALULocalGeometryStorage<GridImp, LocalGeometryObject,2> geoms( myType, false );
        return geoms[ nChild() ];
      }
    }
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: getIndex() const {
    assert( item_ );
    return (*item_).getIndex();
  }

  // forward declararion of struct ElementWrapper
  //template <int codim, int dim, class GridImp>
  //  struct ElementWrapper;

  template<int dim, class GridImp>
  template<int cc>
  inline int ALU2dGridEntity<0,dim,GridImp> :: getSubIndex(int i) const {
    assert(item_ != 0);
    return ElementWrapper<cc, dim, GridImp>::subIndex (grid(), *item_,i);
  }

  //! Provide access to mesh entity i of given codimension. Entities
  //!  are numbered 0 ... count<cc>()-1
  template<int dim, class GridImp>
  template <int cc>
  inline typename ALU2dGridEntity<0,dim,GridImp > ::template Codim<cc> :: EntityPointer
  ALU2dGridEntity<0,dim,GridImp> :: entity (int i) const {
    assert(item_ != 0);
    return ElementWrapper<cc,dim, GridImp>::subEntity (factory(), *item_, i);
  }

  template<int dim, class GridImp>
  template <int cc>
  inline int ALU2dGridEntity<0,dim,GridImp> :: subBoundaryId  ( int i ) const {
    assert(item_ != 0);
    return ElementWrapper<cc, dim, GridImp>::subBoundary (grid(), *item_,i);
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: subIndex(int i, unsigned int codim) const
  {
    assert( item_ != 0 );
    int j = i;
    switch( codim )
    {
    case 0 :
      return ElementWrapper<0, dim, GridImp>::subIndex (grid(), *item_, j);
    case 1 :
      // also apply mapping to generic ref elem by switching edges
      if( item_->numvertices() == 3 )
        j = 2 - i;
      else
        switch (i) { case 0 : j=2;break;
                   case 1 : j=0;break;
                   case 2 : j=3;break;
                   case 3 : j=1;break;}
      // j = ((i^2)>>1) | ((i&1)<<1);
      return ElementWrapper<1, dim, GridImp>::subIndex (grid(), *item_, j);
    case 2 :
      if( item_->numvertices() == 4 )
        switch (i) { case 0 : j=0;break;
                   case 1 : j=1;break;
                   case 2 : j=3;break;
                   case 3 : j=2;break;}
      return ElementWrapper<2, dim, GridImp>::subIndex (grid(), *item_, j);
    default :
      assert( false );
      abort();
    }
    return -1;
  }
  //***************************************************************
  //  Interface for Adaptation
  //***************************************************************
  //! marks an element for refCount refines. if refCount is negative the
  //! element is coarsend -refCount times
  //! mark returns true if element was marked, otherwise false
  template<int dim, class GridImp>
  inline bool ALU2dGridEntity<0,dim,GridImp> :: mark( int refCount ) const
  {
    if( !isLeaf() ) return false;

    // if this assertion is thrown then you try to mark a non leaf entity
    // which is leads to unpredictable results
    assert(item_ != 0);

    // mark for coarsening
    if(refCount < 0)
    {
      if(level() <= 0) return false;
      item_->ALU2DSPACE Refco_el::mark(ALU2DSPACE Refco::crs);
      return true;
    }

    // mark for refinement
    if(refCount > 0)
    {
      item_->ALU2DSPACE Refco_el::mark(ALU2DSPACE Refco::ref);
      return true;
    }

    // mark with none
    item_->ALU2DSPACE Refco_el::mark(ALU2DSPACE Refco::none);
    return true;
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: getMark() const
  {
    assert(item_ != 0);
    if(item_->ALU2DSPACE Refco_el::is(ALU2DSPACE Refco::ref)) return 1;
    if(item_->ALU2DSPACE Refco_el::is(ALU2DSPACE Refco::crs)) return -1;
    assert( item_->ALU2DSPACE Refco_el::is(ALU2DSPACE Refco::none) );
    return 0;
  }

  /*! private methods, but public because of datahandle and template
      arguments of these methods
   */
  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> ::
  setElement(const HElementType &element, int face, int level) const
  {
    item_= const_cast<HElementType *> (&element);
    isLeaf_  = ((*item_).down() == 0);

    geoImpl().unsetUp2Date();
  }

  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> ::
  setElement(const EntitySeed& seed ) const
  {
    setElement( *(seed.item()) ); //, seed.face(), seed.level() );
  }


  //! set actual walk level
  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> :: reset ( int l )
  {
    item_       = 0;
    isLeaf_     = false;

    geoImpl().unsetUp2Date();
  }

  //! set item pointer to NULL
  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> :: removeElement() {
    item_ = 0;
  }

  //! compare 2 entities, which means compare the item pointers
  template<int dim, class GridImp>
  inline bool ALU2dGridEntity<0,dim,GridImp> :: equals ( const ALU2dGridEntity<0,dim,GridImp> & org ) const {
    return (item_ == org.item_);
  }


  //**********************************************************************
  //
  // --EntityPointer
  // --EnPointer
  //
  //**********************************************************************

  //! has to be called when iterator is finished
  template<int cd, class GridImp>
  inline void ALU2dGridEntityPointer<cd, GridImp> :: done()
  {
    // sets entity pointer in the status of an empty entity
    if(entity_)
    {
      entityImp().removeElement();
      factory_.template freeEntity< cd > ( entity_ );
      entity_ = 0;
    }
    seed_.clear();
  }

  template<int cd, class GridImp>
  inline bool ALU2dGridEntityPointer<cd, GridImp> :: equals(const ALU2dGridEntityPointer<cd, GridImp> & i) const
  {
    return seed_ == i.seed_;
  }

  //! update underlying item pointer and set entity
  template<int cd, class GridImp>
  inline void ALU2dGridEntityPointer<cd, GridImp> :: updateEntityPointer(ElementType * item, int face, int level)
  {
    assert(item != 0);
    seed_.set( *item, level, face );

    if( entity_ )
    {
      entityImp().setElement( seed_ );
    }
  }

  //! Constructor for EntityPointer that points to an element
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>::
  ALU2dGridEntityPointer(const FactoryType& factory,
                         const ElementType& item, int face, int level)
    : factory_( factory )
      , seed_( item, level, face )
      , entity_(0)
  { }

  //! Constructor for EntityPointer that points to an element
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>::
  ALU2dGridEntityPointer(const EntityImp& entity)
    : factory_( entity.factory() )
      , seed_( entity.getItem(), entity.level(), entity.getFace() )
      , entity_(0)
  { }

  //! Constructor for EntityPointer that points to an element
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>::
  ALU2dGridEntityPointer(const FactoryType& factory, const EntitySeed& seed)
    : factory_( factory )
      , seed_( seed )
      , entity_(0)
  { }

  //! Constructor for EntityPointer init of Level- and LeafIterator
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>:: ALU2dGridEntityPointer(const FactoryType& factory)
    : factory_( factory )
      , seed_()
      , entity_(0)
  { }

  //! Copy Constructor
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>:: ALU2dGridEntityPointer(const ThisType & org)
    : factory_( org.factory_ )
      , seed_( org.seed_ )
      , entity_(0)
  {  }

  //! Destructor
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>::~ALU2dGridEntityPointer()
  {
    this->done();
  }

  //! dereferencing
  template<int cd, class GridImp>
  inline typename ALU2dGridEntityPointer<cd, GridImp>::Entity &
  ALU2dGridEntityPointer<cd, GridImp>:: dereference() const
  {
    if( ! entity_ )
    {
      entity_ = factory_.template getNewEntity<cd> (level());
      entityImp().setElement( seed_ );
    }
    assert( entity_ );
    return *entity_;
  }

  //! ask for level of entities
  template<int cd, class GridImp>
  inline int ALU2dGridEntityPointer<cd, GridImp>:: level () const
  {
    assert( seed_.level() >= 0 );
    return seed_.level();
  }

  template<int cd, class GridImp>
  inline typename ALU2dGridEntityPointer<cd, GridImp>:: ThisType &
  ALU2dGridEntityPointer<cd, GridImp>:: operator = (const typename ALU2dGridEntityPointer<cd, GridImp>::ThisType & org)
  {
    this->done();
    assert(&factory_ == &org.factory_);
    seed_ = org.seed_; // copy seed
    entity_ = 0; // is set when dereference is called
    return *this;
  }

  template<int cd, class GridImp>
  inline typename ALU2dGridEntityPointer<cd, GridImp>::EntityImp & ALU2dGridEntityPointer<cd, GridImp>::entityImp()
  {
    assert( entity_ );
    return GridImp :: getRealImplementation(*entity_);
  }

  template<int cd, class GridImp>
  inline const typename ALU2dGridEntityPointer<cd, GridImp>:: EntityImp &
  ALU2dGridEntityPointer<cd, GridImp>::entityImp() const {
    assert( entity_ );
    return GridImp :: getRealImplementation(*entity_);
  }

  //********* begin struct ElementWrapper ********************
  //template <int codim, int dim, class GridImp>
  //struct ElementWrapper;
  // partial specialisation for codim
  //
  //--ElementWrapper
  //**********************************************************

  // specialisation for elements
  template<int dim, class GridImp>
  struct ElementWrapper<0,dim, GridImp>
  {
    typedef typename ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType >::HElementType HElementType ;
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    static inline int getElemIndex(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return elem.getIndex();
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return elem.getIndex();
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<0>:: EntityPointer
    subEntity(const FactoryType& factory, const HElementType &elem, int i) {
      //assert(!i);
      return ALU2dGridEntityPointer<0, GridImp > (factory, elem, -1, elem.level());
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return elem.nbbnd(i)->type();

    }
    static inline bool isTheSame(const HElementType * elem, int face, const HElementType * org, int org_face) {
      return (elem == org);
    }
  };

  // specialisation for edges
  template<int dim, class GridImp>
  struct ElementWrapper<1, dim, GridImp>{

    typedef typename ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType >::HElementType HElementType ;
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    static inline int getElemIndex(GridImp & grid, const HElementType &elem, int i)
    {
      assert(i < elem.numvertices() && i >= 0);
      return elem.edge_idx(i);
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i)
    {
      assert(i < elem.numvertices() && i >= 0);
      return elem.edge_idx(i);
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<1>:: EntityPointer
    subEntity(const FactoryType& factory, const HElementType &elem, int i)
    {
      assert(i < elem.numvertices() && i >= 0);
      return ALU2dGridEntityPointer<1, GridImp > (factory, elem, i, elem.level());
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      DUNE_THROW(NotImplemented, "Not yet implemented for this codim!");
      return -1;
    }
    static inline bool isTheSame(const HElementType * elem, int face, const HElementType * org, int org_face)
    {
      if (elem == org)
      {
        if (face == org_face)
          return true;
        else
          return false;
      }
      else {
        if (elem != 0 && org != 0)
          return  (elem->edge_idx(face) == org->edge_idx(org_face));
      }
      return false;
    }
  };

  // specialisation for vertices
  template<int dim, class GridImp>
  struct ElementWrapper<2, dim, GridImp>{

    typedef typename ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType >::HElementType HElementType ;
    typedef typename ALU2dImplInterface< 0, GridImp::dimensionworld, GridImp::elementType >::Type VertexType;
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

    static inline int getElemIndex(GridImp & grid, const VertexType &elem, int) {
      return elem.getIndex();
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      assert(i < elem.numvertices() && i >= 0);
      //return elem.vertex(i)->getIndex();
      return elem.getVertex(i)->getIndex();
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<2>:: EntityPointer
    subEntity(const FactoryType& factory, const HElementType &elem, int i)
    {
      assert(i < elem.numvertices() && i >= 0);
      //return ALU2dGridEntityPointer<2, GridImp > (grid, *(elem.vertex(i)), -1, elem.level());
      return ALU2dGridEntityPointer<2, GridImp > (factory, *(elem.getVertex(i)), -1, elem.level());
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      DUNE_THROW(NotImplemented, "Not yet implemented this codim!");
      return -1;
    }
    static inline bool isTheSame(const VertexType * elem, int face, const VertexType * org, int org_face) {
      return (elem == org);
    }
  };

  //********* end struct ElementWrapper ********************



} //end namespace Dune

#endif
