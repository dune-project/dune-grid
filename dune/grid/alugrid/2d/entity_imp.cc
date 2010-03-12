// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDENTITYIMP_CC
#define DUNE_ALU2DGRIDENTITYIMP_CC

#include "geometry.hh"
#include "grid.hh"
#include <dune/common/exceptions.hh>

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
  ALU2dGridEntity(const GridImp &grid, int level)
    : grid_(grid),
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


  template<>
  inline void ALU2dGridEntity<2,2,ALU2dGrid<2,2> > :: setElement(const HElementType & el, const ALU2DSPACE Vertex & vx)
  {
    item_   = const_cast<ElementType *> (&vx);
    level_  = (*item_).level();

    geoImpl().unsetUp2Date();
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
    : grid_(org.grid_),
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
    return ElementWrapper<cd, dim, GridImp>::getElemIndex (grid_, *item_, face_);
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
  ALU2dGridEntity(const GridImp &grid, int level) :
    grid_(grid)
    , item_(0)
    , geoObj_(GeometryImp())
    , isLeaf_ (false)
  {}

  //! Copy Constructor
  template<int dim, class GridImp>
  inline ALU2dGridEntity<0, dim, GridImp> ::
  ALU2dGridEntity(const ALU2dGridEntity<0,dim,GridImp> & org)
    : grid_(org.grid_),
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
    return EntityPointer(grid_, *(item_->father()));
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
    assert( item_ );
    const int child = this->nChild();
    assert( level() > 0 );

    typedef MakeableInterfaceObject<LocalGeometry> LocalGeometryObject;
    // to be improved, when we using not the refine 8 rule
    if( grid_.nonConform() )
    {
      static ALU2DLocalGeometryStorage<LocalGeometryObject,4> geoms;
      if(!geoms.geomCreated(child))
      {
        typedef typename GridImp::template Codim<0> ::EntityPointer EntityPointer;
        const EntityPointer ep = father();
        geoms.create(grid_, (*ep).geometry(), geometry(), child);
      }
      return geoms[child];
    }
    else
    {
      static ALU2DLocalGeometryStorage<LocalGeometryObject,2> geoms;
      if(!geoms.geomCreated(child))
      {
        typedef typename GridImp::template Codim<0> ::EntityPointer EntityPointer;
        const EntityPointer ep = father();
        geoms.create(grid_,(*ep).geometry(),geometry(),child );
      }
      return geoms[child];
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
    return ElementWrapper<cc, dim, GridImp>::subIndex (grid_, *item_,i);
  }

  //! Provide access to mesh entity i of given codimension. Entities
  //!  are numbered 0 ... count<cc>()-1
  template<int dim, class GridImp>
  template <int cc>
  inline typename ALU2dGridEntity<0,dim,GridImp > ::template Codim<cc> :: EntityPointer
  ALU2dGridEntity<0,dim,GridImp> :: entity (int i) const {
    assert(item_ != 0);
    return ElementWrapper<cc,dim, GridImp>::subEntity (grid_, *item_, i);
  }

  template<int dim, class GridImp>
  template <int cc>
  inline int ALU2dGridEntity<0,dim,GridImp> :: subBoundaryId  ( int i ) const {
    assert(item_ != 0);
    return ElementWrapper<cc, dim, GridImp>::subBoundary (grid_, *item_,i);
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: subIndex(int i, unsigned int codim) const
  {
    assert( item_ != 0 );
    switch( codim )
    {
    case 0 :
      return ElementWrapper<0, dim, GridImp>::subIndex (grid_, *item_, i);
    case 1 :
      // also apply mapping to generic ref elem by switching edges
      return ElementWrapper<1, dim, GridImp>::subIndex (grid_, *item_, 2-i);
    case 2 :
      return ElementWrapper<2, dim, GridImp>::subIndex (grid_, *item_, i);
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
      item_->Refco_el::mark(ALU2DSPACE Refco::crs);
      return true;
    }

    // mark for refinement
    if(refCount > 0)
    {
      item_->Refco_el::mark(ALU2DSPACE Refco::ref);
      return true;
    }

    // mark with none
    item_->Refco_el::mark(ALU2DSPACE Refco::none);
    return true;
  }

  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: getMark() const
  {
    assert(item_ != 0);
    if(item_->Refco_el::is(ALU2DSPACE Refco::ref)) return 1;
    if(item_->Refco_el::is(ALU2DSPACE Refco::crs)) return -1;
    assert( item_->Refco_el::is(ALU2DSPACE Refco::none) );
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
    item_ = 0;
    face_ = -1; // set face to non-valid value
    compactify();
  }

  template<int cd, class GridImp>
  inline void ALU2dGridEntityPointer<cd, GridImp> :: compactify()
  {
    // sets entity pointer in the status of an empty entity
    if(entity_)
    {
      entityImp().removeElement();
      grid_.freeEntity( entity_ );
      entity_ = 0;
    }
  }

  template<int cd, class GridImp>
  inline bool ALU2dGridEntityPointer<cd, GridImp> :: equals(const ALU2dGridEntityPointer<cd, GridImp> & i) const
  {
    return ElementWrapper<cd,dim, GridImp>::isTheSame (item_, face_, i.item_, i.face_);
  }

  //! update underlying item pointer and set entity
  template<int cd, class GridImp>
  inline void ALU2dGridEntityPointer<cd, GridImp> :: updateEntityPointer(ElementType * item, int face, int level) {
    assert(item != 0);
    item_ = item;
    assert(item_);

    face_= face;
    level_ = level;
    if( entity_ )
    {
      entityImp().setElement( *item_, face_, level_);
    }
  }

  //! update underlying item pointer and set entity
  //! specialization for codim 0
  template<>
  inline void ALU2dGridEntityPointer<0, const ALU2dGrid<2,2> > :: updateEntityPointer(ElementType * item, int , int ) {
    assert(item != 0);
    item_ = item;
    assert(item_);

    if( entity_ )
    {
      entityImp().setElement( *item_, -1 , -1 );
    }
  }


  //! Constructor for EntityPointer that points to an element
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>:: ALU2dGridEntityPointer(const GridImp & grid,
                                                                      const ElementType & item, int face, int level)
    : grid_(grid)
      , item_(const_cast<ElementType *>(&item))
      , level_(level)
      , face_(face)
      , entity_(0)
  { }

  //! Constructor for EntityPointer that points to an element
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>::
  ALU2dGridEntityPointer(const EntityImp& entity)
    : grid_(entity.grid())
      , item_(& entity.getItem() )
      , level_(entity.level())
      , face_(entity.getFace())
      , entity_(0)
  { }

  //! Constructor for EntityPointer init of Level- and LeafIterator
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>:: ALU2dGridEntityPointer(const GridImp & grid)
    : grid_(grid)
      , item_(0)
      , level_(-1)
      , face_(-1)
      , entity_(0)
  { }

  //! Copy Constructor
  template<int cd, class GridImp>
  inline ALU2dGridEntityPointer<cd, GridImp>:: ALU2dGridEntityPointer(const ThisType & org)
    : grid_(org.grid_)
      , item_(org.item_)
      , level_(org.level_)
      , face_(org.face_)
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
    assert( item_ );
    if( !entity_ )
    {
      entity_ = grid_.template getNewEntity<cd> (level());
      entityImp().setElement(*item_, face_, level());
    }
    assert( entity_ );
    return *entity_;
  }

  //! ask for level of entities
  template<int cd, class GridImp>
  inline int ALU2dGridEntityPointer<cd, GridImp>:: level () const
  {
    assert( item_ );
    if (level_ == -1)
    {
      if (cd == 2)
      {
        // ????
        level_ = item_->level()+1;
      }
      else
      {
        level_ = item_->level();
      }
    }
    return level_;
  }

  //! ask for level of entities
  template<>
  inline int ALU2dGridEntityPointer<0, const ALU2dGrid<2,2> >:: level () const
  {
    assert( item_ );
    return item_->level();
  }

  template<int cd, class GridImp>
  inline typename ALU2dGridEntityPointer<cd, GridImp>:: ThisType &
  ALU2dGridEntityPointer<cd, GridImp>:: operator = (const typename ALU2dGridEntityPointer<cd, GridImp>::ThisType & org)
  {
    this->done();
    entity_ = 0;
    assert(&grid_ == &org.grid_);
    item_  = org.item_;
    face_  = org.face_;
    level_ = org.level_;
    return *this;
  }

  template<int cd, class GridImp>
  inline typename ALU2dGridEntityPointer<cd, GridImp>::EntityImp & ALU2dGridEntityPointer<cd, GridImp>::entityImp()
  {
    assert( entity_ );
    return grid_.getRealImplementation(*entity_);
  }

  template<int cd, class GridImp>
  inline const typename ALU2dGridEntityPointer<cd, GridImp>:: EntityImp &
  ALU2dGridEntityPointer<cd, GridImp>::entityImp() const {
    assert( entity_ );
    return grid_.getRealImplementation(*entity_);
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
  struct ElementWrapper<0,dim, GridImp>{

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    static inline int getElemIndex(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return elem.getIndex();
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return elem.getIndex();
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<0>:: EntityPointer
    subEntity(GridImp & grid, const HElementType &elem, int i) {
      //assert(!i);
      return ALU2dGridEntityPointer<0, GridImp > (grid, elem, -1, elem.level());
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

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    static inline int getElemIndex(GridImp & grid, const HElementType &elem, int i)
    {
      assert(i < 3 && i >= 0);
      return elem.edge_idx(i);
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i)
    {
      assert(i < 3 && i >= 0);
      return elem.edge_idx(i);
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<1>:: EntityPointer
    subEntity(GridImp & grid, const HElementType &elem, int i)
    {
      assert(i < 3 && i >= 0);
      return ALU2dGridEntityPointer<1, GridImp > (grid, elem, i, elem.level());
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

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    static inline int getElemIndex(GridImp & grid, const ALU2DSPACE Vertex &elem, int) {
      return elem.getIndex();
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      //return elem.vertex(i)->getIndex();
      return elem.getVertex(i)->getIndex();
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<2>:: EntityPointer
    subEntity(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      //return ALU2dGridEntityPointer<2, GridImp > (grid, *(elem.vertex(i)), -1, elem.level());
      return ALU2dGridEntityPointer<2, GridImp > (grid, *(elem.getVertex(i)), -1, elem.level());
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      DUNE_THROW(NotImplemented, "Not yet implemented this codim!");
      return -1;
    }
    static inline bool isTheSame(const ALU2DSPACE Vertex * elem, int face, const ALU2DSPACE Vertex * org, int org_face) {
      return (elem == org);
    }
  };

  //********* end struct ElementWrapper ********************



} //end namespace Dune

#endif
