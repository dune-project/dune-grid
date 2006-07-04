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
      father_(0),
      geoObj_(GeometryImp()),
      geoImp_(grid_.getRealImplementation(geoObj_)),
      builtgeometry_(false),
      level_(level),
      face_(-1),
      localFCoordCalced_(false) {}


  template<int cd, int dim, class GridImp>
  inline void ALU2dGridEntity<cd,dim,GridImp>:: setElement(const ElementType &element, int face, int level) const {
    item_= const_cast<ElementType *> (&element);
    builtgeometry_=false;
    //level_ = item_->level();
    level_ = level;
    face_ = face;
    localFCoordCalced_=false;
  }


  template<>
  inline void ALU2dGridEntity<2,2,ALU2dGrid<2,2> > :: setElement(const HElementType & el, const ALU2DSPACE Vertex & vx) {
    item_   = const_cast<ElementType *> (&vx);
    level_  = (*item_).level();
    father_ = (&el);
    builtgeometry_=false;
    localFCoordCalced_=false;
  }

  //! set item pointer to NULL
  template<int cd, int dim, class GridImp>
  inline void ALU2dGridEntity<cd,dim,GridImp> :: removeElement() {
    item_ = 0;
    father_=0;
  }

  //! Copy Constructor
  template<int cd, int dim, class GridImp>
  inline ALU2dGridEntity<cd, dim, GridImp> ::
  ALU2dGridEntity(const ALU2dGridEntity<cd,dim,GridImp> & org)
    : grid_(org.grid_),
      item_(org.item_),
      father_(org.father_),
      geoObj_(GeometryImp()),
      geoImp_(grid_.getRealImplementation(geoObj_)),
      builtgeometry_(false),
      level_(org.level_),
      face_(org.face_),
      localFCoordCalced_(org.localFCoordCalced_),
      localFatherCoords_() { }

  //! geometry of this entity
  template<int cd, int dim, class GridImp>
  inline const typename ALU2dGridEntity<cd, dim, GridImp> :: Geometry & ALU2dGridEntity<cd, dim, GridImp> ::
  geometry () const {
    if(!builtgeometry_) builtgeometry_ = geoImp_.builtGeom(*item_,face_);
    assert(builtgeometry_ == true);
    return geoObj_;
  }

  template<int cd, int dim, class GridImp>
  inline int ALU2dGridEntity<cd,dim,GridImp > :: getIndex() const {
    assert(item_ != 0);
    return ElementWrapper<cd, dim, GridImp>::getElemIndex (grid_, *item_, face_);
    //     if (cd != 1)
    //       return item_->getIndex();
    //     else
    //       return item_->edge_idx(face_);
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


  template<int cd, int dim, class GridImp>
  inline typename ALU2dGridEntity<cd,dim,GridImp>::EntityPointer
  ALU2dGridEntity<cd,dim,GridImp>:: ownersFather() const {
    assert(cd == dim); // this method only exists for codim == dim
    if( !father_ )
    {
      dwarn << "No Father for given Entity! \n";
      return ALU2dGridEntityPointer<0,GridImp> (grid_,(*father_));
    }
    return ALU2dGridEntityPointer<0,GridImp> (grid_,(*father_));
  }

  template<int cd, int dim, class GridImp>
  inline FieldVector<alu2d_ctype, dim> &
  ALU2dGridEntity<cd,dim,GridImp>:: positionInOwnersFather() const {
    assert( cd == dim );
    if(!localFCoordCalced_)
    {
      EntityPointer vati = this->ownersFather();
      localFatherCoords_ = (*vati).geometry().local( this->geometry()[0] );
      localFCoordCalced_ = true;
    }
    return localFatherCoords_;
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
    , geoImp_(grid_.getRealImplementation(geoObj_))
    , builtgeometry_(false)
    //, index_()
    , walkLevel_ (level)
    , isLeaf_ (false) {}

  //! Copy Constructor
  template<int dim, class GridImp>
  inline ALU2dGridEntity<0, dim, GridImp> ::
  ALU2dGridEntity(const ALU2dGridEntity<0,dim,GridImp> & org)
    : grid_(org.grid_),
      item_(org.item_),
      geoObj_(GeometryImp()),
      geoImp_(grid_.getRealImplementation(geoObj_)),
      builtgeometry_(false),
      walkLevel_(org.walkLevel_),
      isLeaf_(org.isLeaf_) { }

  //! level of this element
  template<int dim, class GridImp>
  inline int ALU2dGridEntity<0,dim,GridImp> :: level () const {
    assert( item_ );
    return (*item_).level();
  }

  //! geometry of this entity
  template<int dim, class GridImp>
  inline const typename ALU2dGridEntity<0, dim, GridImp> :: Geometry & ALU2dGridEntity<0,dim,GridImp> :: geometry () const {
    assert(item_ != 0);
    if(!builtgeometry_) builtgeometry_ = geoImp_.builtGeom(*item_,-1);

    assert(builtgeometry_ == true);
    return geoObj_;
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
  inline int ALU2dGridEntity<0, dim, GridImp> :: nChild() const {
    return item_->childNr();
  }

  // singletons of geometry in father geometries
  // GeometryType schould be of type Dune::Geometry
  /*
     template <class GeometryType>
     static inline GeometryType &
     getGeometryInFather(const int child, const int orientation = 1)
     {
     typedef typename GeometryType :: ImplementationType GeometryImp;
     static GeometryType child0       (GeometryImp(0,1)); // child 0
     static GeometryType child1_plus  (GeometryImp(1,1)); // child 1
     static GeometryType child1_minus (GeometryImp(1,-1)); // child 1, orientation < 0

     if(child == 0) return child0;
     if(child == 1) return (orientation > 0) ? child1_plus : child1_minus;

     DUNE_THROW(NotImplemented,"wrong number of child given!");
     return child0;
     }
   */

  template<int dim, class GridImp>
  inline const typename ALU2dGridEntity<0, dim, GridImp>:: Geometry & ALU2dGridEntity<0,dim,GridImp> ::
  geometryInFather () const
  {
    assert( item_ );
    const int child = this->nChild();
    typedef typename Geometry::ImplementationType GeometryImp;
    // to be improved, when we using not the refine 8 rule
    static ALU2DLocalGeometryStorage<Geometry,4> geoms;
    if(!geoms.geomCreated(child))
    {
      typedef typename GridImp::template Codim<0> ::EntityPointer EntityPointer;
      const EntityPointer ep = father();
      geoms.create(grid_,(*ep).geometry(),geometry(),child );
    }
    return geoms[child];

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
    return ElementWrapper<cc,dim, GridImp>::subEntity (grid_, *item_,i);
  }

  template<int dim, class GridImp>
  template <int cc>
  inline int ALU2dGridEntity<0,dim,GridImp> :: subBoundaryId  ( int i ) const {
    assert(item_ != 0);
    return ElementWrapper<cc, dim, GridImp>::subBoundary (grid_, *item_,i);
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
    //assert(item_ != 0);

    // if this assertion is thrown then you try to mark a non leaf entity
    // which is leads to unpredictable results
    assert( isLeaf() );

    // mark for coarsening
    if(refCount < 0) {
      if(level() <= 0) return false;
      item_->Refco_el::mark(ALU2DSPACE Refco::crs);
    }
    // mark for refinement
    if(refCount > 0) {
      item_->Refco_el::mark(ALU2DSPACE Refco::ref);
    }
    else
      return false;
    return true;
  }

  /*! private methods, but public because of datahandle and template
      arguments of these methods
   */
  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> :: setElement(const HElementType &element, int face, int level) const {
    item_= const_cast<HElementType *> (&element);
    builtgeometry_=false;
    //index_   = -1;
    //level_   = (*item_).level();
    //glIndex_ = (*item_).getIndex();
    isLeaf_  = ((*item_).down() == 0);
  }

  //! set actual walk level
  template<int dim, class GridImp>
  inline void ALU2dGridEntity<0,dim,GridImp> :: reset ( int l ){
    assert( walkLevel_ >= 0 );

    item_       = 0;
    builtgeometry_ = false;
    walkLevel_     = l;
    //glIndex_    = -1;
    //level_      = -1;
    isLeaf_     = false;
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
  inline void ALU2dGridEntityPointer<cd, GridImp> :: done() {
    item_ = 0;
    face_ = -1; // set face to non-valid value
    // sets entity pointer in the status of an empty entity
    if(entity_)
    {
      entityImp().removeElement();
      //grid_.freeEntity( entity_ );
      delete entity_;
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
    if( item_ && entity_ ) {
      entityImp().setElement( *item_, face_, level_);
    }
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
      return ALU2dGridEntityPointer<0, GridImp > (grid, elem);
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

    static inline int getElemIndex(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      return elem.edge_idx(i);
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      return elem.edge_idx(i);
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<1>:: EntityPointer
    subEntity(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      return ALU2dGridEntityPointer<1, GridImp > (grid, elem, i);
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      DUNE_THROW(NotImplemented, "Not yet implemented for this codim!");
    }
    static inline bool isTheSame(const HElementType * elem, int face, const HElementType * org, int org_face) {
      if (elem == org) {
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

    static inline int getElemIndex(GridImp & grid, const ALU2DSPACE Vertex &elem, int i) {
      assert(i < 3 && i >= 0);
      //cout << "  getElemIndex: " << elem.getIndex() << endl;
      //cout << "  coord: " << elem.coord()[0] << " " << elem.coord()[1] << endl;
      return elem.getIndex();
    }
    static inline int subIndex(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      //cout << "  subIndex: " << elem.vertex(i)->getIndex() << endl;
      //cout << "  coord: " << elem.vertex(i)->coord()[0] << " " << elem.vertex(i)->coord()[1] << endl;
      return elem.vertex(i)->getIndex();
    }
    static inline typename ALU2dGridEntity<0,dim,GridImp > :: template Codim<2>:: EntityPointer
    subEntity(GridImp & grid, const HElementType &elem, int i) {
      assert(i < 3 && i >= 0);
      return ALU2dGridEntityPointer<2, GridImp > (grid, *(elem.vertex(i)));
    }
    static inline int subBoundary(GridImp & grid, const HElementType &elem, int i) {
      DUNE_THROW(NotImplemented, "Not yet implemented this codim!");
    }
    static inline bool isTheSame(const ALU2DSPACE Vertex * elem, int face, const ALU2DSPACE Vertex * org, int org_face) {
      return (elem == org);
    }
  };

  //********* end struct ElementWrapper ********************



} //end namespace Dune

#endif
