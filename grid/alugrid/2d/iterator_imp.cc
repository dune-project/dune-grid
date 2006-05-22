// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "geometry.hh"
#include "entity.hh"
#include "grid.hh"
//#include "faceutility.hh"


#ifndef DUNE_ALU2DGRID_ITERATOR_IMP_CC
#define DUNE_ALU2DGRID_ITERATOR_IMP_CC


namespace Dune {
  //********************************************************************
  //
  //  --ALU2dGridIntersectionIterator
  //  --IntersectionIterator
  //
  //********************************************************************


  // --IntersectionIterator
  //! Constructor
  template<class GridImp>
  inline ALU2dGridIntersectionIterator<GridImp> ::
  ALU2dGridIntersectionIterator(const GridImp & grid, const HElementType* el, int wLevel, bool end) :
    //geoProvider_(connector_),
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(grid),
    item_(0),
    neigh_(0),
    nFaces_(3),
    walkLevel_(wLevel),
    index_(0),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    done_(end)
  {
    if (!end)
    {
      assert(walkLevel_ >= 0);
      setFirstItem(*el,wLevel);
    }
    else
    {
      done();
    }
  }

  template<class GridImp>
  inline ALU2dGridIntersectionIterator<GridImp> ::
  ALU2dGridIntersectionIterator(const GridImp & grid, int wLevel) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(grid),
    item_(0),
    neigh_(0),
    nFaces_(3),
    walkLevel_(wLevel),
    index_(0),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    done_(true) ,
    calledOnLeaf_(false),
    visitedLeaf_(false)
  {}


  //! The copy constructor
  template<class GridImp>
  inline ALU2dGridIntersectionIterator<GridImp> ::
  ALU2dGridIntersectionIterator(const ALU2dGridIntersectionIterator<GridImp> & org) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(org.grid_),
    item_(org.item_),
    neigh_(org.neigh_),
    nFaces_(org.nFaces_),
    walkLevel_(org.walkLevel_),
    index_(org.index_),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    //outerNormal_(org.unitOuterNormal_), call default constructor of these
    //unitOuterNormal_(org.unitOuterNormal_),
    done_(org.done_),
    calledOnLeaf_(org.calledOnLeaf_),
    visitedLeaf_(org.visitedLeaf_)
  {}

  //! The copy constructor
  template<class GridImp>
  inline void
  ALU2dGridIntersectionIterator<GridImp> ::
  assign(const ALU2dGridIntersectionIterator<GridImp> & org)
  {
    assert( &grid_ == &org.grid_);
    item_      = org.item_;
    neigh_     = org.neigh_;
    index_     = org.index_;
    nFaces_    = org.nFaces_;
    walkLevel_ = org.walkLevel_;
    generatedGlobalGeometry_ = false;
    generatedLocalGeometries_ = false;
    done_ = org.done_;
    calledOnLeaf_ = org.calledOnLeaf_;
    visitedLeaf_ = org.visitedLeaf_;
  }

  //! check whether entities are the same or whether iterator is done
  template<class GridImp>
  inline bool ALU2dGridIntersectionIterator<GridImp> ::
  equals (const ALU2dGridIntersectionIterator<GridImp> & i) const
  {
    return ((item_ == i.item_) &&
            (done_ == i.done_));

  }


  //! increment iterator
  template<class GridImp>
  inline void ALU2dGridIntersectionIterator<GridImp> :: increment ()
  {
    if (index_ >= nFaces_) {
      done();
      return ;
    }
    else {
      if(!calledOnLeaf_) // we don't want leaf intersections
      {
        ++index_;
        if (index_ >= nFaces_)
        {
          done();
          return ;
        }
        neigh_ = item_->nbel(index_); //! maybe NULL, check boundary()
        if (neigh_ == 0) // on boundary do nothing
          return;
        else {
          int nlevel = neigh_->level();
          if (nlevel == walkLevel_) // intersection is on desired level
            return;
          if (nlevel < walkLevel_) { // then we don't have an intersection on walkLevel_
            increment();
            return;
          }
          else {
            while (neigh_ != 0 && nlevel > walkLevel_) {
              neigh_ = neigh_->father(); // maybe NULL?
              nlevel = neigh_->level();
            }
            return;
          }
        }
      }
      else { // now return leaf elements too
        if (!visitedLeaf_)
        {
          ++index_;
          if (index_ >= nFaces_)
          {
            done();
            return ;
          }
          neigh_ = item_->nbel(index_); //! maybe NULL, check boundary()
          visitedLeaf_ = true;
        }
        else {
          if (neigh_ == 0) { // on boundary do nothing
            visitedLeaf_ = false;
            return;
          }
          else {
            int nlevel = neigh_->level();
            if (nlevel == walkLevel_) { // intersection is on desired level
              visitedLeaf_ = false;
              return;
            }
            if (nlevel < walkLevel_) { // then we don't have an intersection on walkLevel_
              visitedLeaf_ = false;
              increment();
              return;
            }
            else {
              assert(nlevel > 0);
              while (neigh_ != 0 && nlevel > walkLevel_) {
                neigh_ = neigh_->father(); // maybe NULL?
                assert(neigh_ != 0);
                nlevel = neigh_->level();
                assert(nlevel > 0);
              }
              assert(neigh_ != 0);
              visitedLeaf_ = false;
              return;
            }
          }
        }
      }
    }
  }

  /*
        if (index_ >= nFaces_) {
        done();
        return ;
      }
     ++index_;
      if (index_ >= nFaces_)
      {
        done();
        return ;
      }
      neigh_ = item_->nbel(index_);  //! maybe NULL
     }
   */


  //! return level of inside() entitiy
  template<class GridImp>
  inline int ALU2dGridIntersectionIterator<GridImp> :: level () const {
    return (*item_).level();
  }


  //! set interator to end iterator
  template<class GridImp>
  inline void ALU2dGridIntersectionIterator<GridImp> :: done () {
    done_ = true;
    item_ = 0;
    neigh_ = 0;
    index_= nFaces_;
  }


  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  template<class EntityType>
  inline void ALU2dGridIntersectionIterator<GridImp> ::
  first(const EntityType & en, int wLevel)
  {
    setFirstItem(en.getItem(),wLevel);
  }

  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  inline void ALU2dGridIntersectionIterator<GridImp> :: setFirstItem(const HElementType & elem, int wLevel) {

    item_      = const_cast<HElementType *> (&elem);
    assert( item_ );
    walkLevel_ = wLevel;
    done_ = false;
    calledOnLeaf_= (elem.down() == 0);
    visitedLeaf_ = false;

    index_ = 0;
    //while (!item_->nbel(index_) && index_ < 3)
    //  ++index_;
    neigh_ = item_->nbel(index_); //! maybe NULL
  }


  //! return true if intersection is with boundary
  template<class GridImp>
  inline bool ALU2dGridIntersectionIterator<GridImp> :: boundary() const {
    return (neigh_ == 0);
    //ALU2DSPACE Thinelement * nb = item_->neighbour(index_);
    //return nb->thinis(ALU2DSPACE Thinelement::bndel_like);
  }

  template<class GridImp>
  inline int ALU2dGridIntersectionIterator<GridImp> :: boundaryId() const
  {
    int isBoundary=0;
    assert(item_);
    if(item_->nbbnd(index_) != 0)
      isBoundary = item_->nbbnd(index_)->type();
    return isBoundary;
  }


  //! return true if intersection is with neighbor on this level
  template<class GridImp>
  inline bool ALU2dGridIntersectionIterator<GridImp> :: neighbor () const {
    return !(this->boundary());
  }

  //! return true if intersection is with neighbor on this level
  template<class GridImp>
  inline bool ALU2dGridIntersectionIterator<GridImp> :: leafNeighbor () const {
    if (!boundary() && calledOnLeaf_)
      return (neigh_->down()==0);
    else
      return false;
  }

  //! return true if intersection is with neighbor on this level
  template<class GridImp>
  inline bool ALU2dGridIntersectionIterator<GridImp> :: levelNeighbor () const {
    if (!boundary() && neigh_ != 0)
      return (neigh_->level()==walkLevel_);
    else
      return false;
  }

  //! return EntityPointer to the Entity on the inside of this intersection.
  template<class GridImp>
  inline typename ALU2dGridIntersectionIterator<GridImp> :: EntityPointer
  ALU2dGridIntersectionIterator<GridImp> :: inside() const {
    assert( item_ );
    return EntityPointer(grid_, *item_);
  }


  //! return EntityPointer to the Entity on the outside of this intersection.
  template<class GridImp>
  inline typename ALU2dGridIntersectionIterator<GridImp> :: EntityPointer
  ALU2dGridIntersectionIterator<GridImp> :: outside() const {
    assert(!boundary());
    assert( neigh_ );
    return EntityPointer(grid_, *neigh_);
  }


  //! local number of codim 1 entity in self where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionIterator<GridImp> :: numberInSelf () const {
    return index_;
  }

  //! local number of codim 1 entity in neighbor where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionIterator<GridImp> :: numberInNeighbor () const {
    if (neigh_->down() == 0)
      return item_->opposite(index_);
    else {
      int i=0;
      assert(index_ < 3);
      assert(neigh_->level() == walkLevel_);
      //for (int i=0; i < 3; ++i)
      //  std::cout << index << "  " <<  item_->vertex(i)->coord()[0] << "  " << item_->vertex(i)->coord()[1] << endl ;
      //for (int i=0; i < 3; ++i)
      //  std::cout << neigh_->vertex(i)->coord()[0] << "  " << neigh_->vertex(i)->coord()[1] << endl;


      EntityPointer item(grid_, *item_);
      EntityPointer neigh(grid_, *neigh_);
      while(i<3 && grid_.levelIndexSet(walkLevel_).template subIndex<1>(item.dereference(), index_)!= grid_.levelIndexSet(walkLevel_).template subIndex<1>(neigh.dereference(),i)) {
        ++i;
      }
      //assert(i<3);
      return i;
    }

  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionIterator<GridImp>::NormalType &
  ALU2dGridIntersectionIterator<GridImp> :: outerNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    assert(item_ != 0);

    double dummy[2];
    item_->outernormal(index_, dummy);
    outerNormal_[0] = dummy[0];
    outerNormal_[1] = dummy[1];

    return outerNormal_;
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionIterator<GridImp>::NormalType &
  ALU2dGridIntersectionIterator<GridImp> :: integrationOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    return this->outerNormal(local);
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionIterator<GridImp>::NormalType &
  ALU2dGridIntersectionIterator<GridImp> :: unitOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    unitOuterNormal_ = this->outerNormal(local);
    unitOuterNormal_ *= (1.0/unitOuterNormal_.two_norm());
    return unitOuterNormal_;
  }


  template<class GridImp>
  inline const typename ALU2dGridIntersectionIterator<GridImp>::LocalGeometry&
  ALU2dGridIntersectionIterator<GridImp> ::intersectionSelfLocal () const {
    assert(item_ != 0);
    typename ALU2dGridIntersectionIterator<GridImp> :: EntityPointer ep = inside();
    this->grid_.getRealImplementation(intersectionSelfLocal_).builtLocalGeom( ep.dereference().geometry() , intersectionGlobal() );
    //, intersectionSelfLocal_, item_, index_);
    return intersectionSelfLocal_;
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionIterator<GridImp>::LocalGeometry&
  ALU2dGridIntersectionIterator<GridImp> :: intersectionNeighborLocal () const {
    assert(item_ != 0 && neigh_ != 0);
    typename ALU2dGridIntersectionIterator<GridImp> :: EntityPointer ep = outside();
    this->grid_.getRealImplementation(intersectionNeighborLocal_).builtLocalGeom( ep.dereference().geometry() , intersectionGlobal()); //, intersectionNeighborLocal_, neigh_ , numberInNeighbor());
    return intersectionNeighborLocal_;
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionIterator<GridImp>::Geometry&
  ALU2dGridIntersectionIterator<GridImp> ::intersectionGlobal () const {
    assert(item_ != 0);
    this->grid_.getRealImplementation(intersectionGlobal_).builtGeom(*item_, index_);
    return intersectionGlobal_;
  }


  //********* begin struct CheckElement ********************
  //template<int cc, PartitionIteratorType, class GridImp>
  //struct CheckElementType;
  //
  //********************************************************

  template<int cc, PartitionIteratorType, class GridImp>
  struct CheckElementType;

  // specialisation for elements
  template<PartitionIteratorType pitype, class GridImp>
  struct CheckElementType<0,pitype,GridImp>{
    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    static inline int checkFace(HElementType & item, int & face, int level, bool isLeafIterator) {
      return 1;
    }
  };

  // specialisation for edges
  template<PartitionIteratorType pitype, class GridImp>
  struct CheckElementType<1,pitype,GridImp>{
    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    //static inline int checkFace(typename ALU2dGridLeafIterator<1,pitype,GridImp>::ElementType & item, int & face) {
    static inline int checkFace(HElementType & item, int & face, int level, bool isLeafIterator) {
      assert(face>=0);
      assert(isLeafIterator);

      while (face < 3) {
        if(item.normaldir(face)==1) {
          return 0;
        }
        else ;
        ++face;
      }
      return 1;
    }
  };

  // specialisation for vertices
  template<PartitionIteratorType pitype, class GridImp>
  struct CheckElementType<2,pitype,GridImp>{
    static inline int checkFace(ALU2DSPACE Vertex & item, int & face, int level, bool isLeafIterator) {
      return 1;
    }
  };
  //********* end struct CheckElement ********************


  //********************************************************************
  //
  //  --TreeIterator
  //
  //********************************************************************

  //! constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline TreeIterator<cdim, pitype, GridImp> ::
  TreeIterator(const GridImp & grid, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(-1),
    face_(0),
    isCopy_(0),
    elem_(0),
    iter_(),
    isLeafIterator_(true)
  {
    if(!end) {
      IteratorType * it = new IteratorType(const_cast<ALU2DSPACE Hmesh&>(grid.myGrid()));
      iter_.store( it );

      IteratorType & iter = *it;

      iter->first();
      if((!iter->done()))
      {
        elem_ = &(iter->getitem());
        this->updateEntityPointer(elem_, face_, elem_->level());
        if(cdim==1)
          increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline TreeIterator<cdim, pitype, GridImp> ::
  TreeIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    face_(0),
    isCopy_(0),
    elem_(0),
    iter_(),
    isLeafIterator_(false)
  {

    if(!end) {
      IteratorType * it = new IteratorType(const_cast<ALU2DSPACE Hmesh&>(grid.myGrid()), level_);

      iter_.store( it );

      IteratorType & iter = *it;

      iter->first();
      if((!iter->done()))
      {
        elem_ = &(iter->getitem());
        this->updateEntityPointer(elem_, face_, level_);
        if(cdim==1)
          increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }


  //! copy constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline TreeIterator<cdim, pitype, GridImp> ::
  TreeIterator(const TreeIterator<cdim,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , face_(org.face_)
      , isCopy_(org.isCopy_+1)
      , elem_(org.elem_)
      , iter_ ( org.iter_ )
      , isLeafIterator_(org.isLeafIterator_)
  {
    // don't copy a copy of a copy of a copy of a copy
    assert( org.isCopy_ < 3 );
  }

  //! prefix increment
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline void TreeIterator<cdim, pitype, GridImp> :: increment () {
    //assert(face_>=0);
    if(endIter_)
      return ;

    IteratorType & iter = iter_.operator *();
    int goNext = CheckElementType<cdim,pitype,GridImp>::checkFace(*(this->item_), face_, level_, isLeafIterator_);

    if (goNext) {
      if (cdim ==1) {
        assert(face_==3);
        iter->next();
        if(iter->done()) {
          endIter_ = true;
          face_= 0;
          this->done();
          return ;
        }
        face_=0;
        elem_ = &(iter->getitem());
        this->updateEntityPointer(elem_, face_);
        increment();
        return;
      }
      else {
        iter->next();
        face_=0;
      }
    }

    if(!goNext || cdim!= 1) {
      if(iter->done()) {
        endIter_ = true;
        face_= 0;
        this->done();
        return ;
      }

      elem_ = &(iter->getitem());
      this->updateEntityPointer(elem_, face_);
      ++face_;
    }
  }



  //********************************************************************
  //
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************


  //! Constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU2dGridLeafIterator(const GridImp & grid, bool end) :
    TreeIterator<cdim,pitype,GridImp> (grid, end) {  }


  //! copy constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU2dGridLeafIterator(const ALU2dGridLeafIterator<cdim, pitype, GridImp> & org) :
    TreeIterator<cdim,pitype,GridImp> (org) {  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator
  //
  //********************************************************************


  //! constructor
  template<int cd, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<cd, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    TreeIterator<cd,pitype,GridImp> (grid, level, end) {  }

  //! copy constructor
  template<int cd, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<cd, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<cd,pitype,GridImp> & org)
    :   TreeIterator<cd,pitype,GridImp> (org) {  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator, specialized for cd=2
  //
  //********************************************************************

  //! constructor for cd=2
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    face_(0),
    isCopy_(0),
    nrOfVertices_(grid.size(2)),
    iter_() {

    indexList = new double[nrOfVertices_];
    for (int i = 0; i < nrOfVertices_; ++i)
      indexList[i]= 0;

    if(!end) {
      IteratorType * it = new IteratorType(const_cast<ALU2DSPACE Hmesh&>(grid.myGrid()), level_);
      iter_.store( it );
      IteratorType & iter = *it;
      iter->first();
      if((!iter->done()))
      {
        item_ = &iter->getitem();
        vertex_ = item_->vertex(face_);
        indexList[vertex_->getIndex()] = 1;
        this->updateEntityPointer(vertex_, face_, level_);
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor for cd=2
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<2,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , face_(org.face_)
      , isCopy_(org.isCopy_+1)
      , nrOfVertices_(org.nrOfVertices_)
      , item_(org.item_)
      , vertex_(org.vertex_)
      , iter_ ( org.iter_ )
  {

    indexList = new double[nrOfVertices_];
    for (int i = 0; i < nrOfVertices_; ++i)
      indexList[i] = org.indexList[i];

    // don't copy a copy of a copy of a copy of a copy
    assert( org.isCopy_ < 3 );
  }

  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLevelIterator<2, pitype, GridImp> :: increment () {

    if(endIter_)
      return ;

    IteratorType & iter = iter_.operator *();

    assert(face_>=0);
    int goNext = 1;
    item_ = &iter->getitem();
    while (face_ < 3) {
      vertex_ = item_->vertex(face_);
      int idx = vertex_->getIndex();
      if(!indexList[idx]) {
        indexList[idx]=1;
        goNext = 0;
        break;
      }
      ++face_;
    }

    if (goNext) {
      assert(face_==3);
      iter->next();
      if(iter->done()) {
        endIter_ = true;
        face_= 0;
        this->done();
        return ;
      }
      face_=0;
      item_ = &iter->getitem();
      vertex_ = item_->vertex(face_);
      this->updateEntityPointer(vertex_, face_, level_);
      increment();
      return;
    }

    if(iter->done()) {
      endIter_ = true;
      face_= 0;
      this->done();
      return ;
    }
    item_ = &iter->getitem();
    vertex_ = item_->vertex(face_);
    this->updateEntityPointer(vertex_, face_, level_);
    ++face_;
  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator, specialized for cd=1
  //
  //********************************************************************

  //! constructor for cd=1
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    face_(0),
    isCopy_(0),
    nrOfEdges_(grid.hierSetSize(1)),
    iter_() {

    indexList = new double[nrOfEdges_];
    for (int i = 0; i < nrOfEdges_; ++i)
      indexList[i]= 0;

    if(!end) {
      IteratorType * it = new IteratorType(const_cast<ALU2DSPACE Hmesh&>(grid.myGrid()), level_);
      iter_.store( it );
      IteratorType & iter = *it;
      iter->first();
      if((!iter->done()))
      {
        elem_ = &(iter->getitem());
        this->updateEntityPointer(elem_, face_, level_);
        increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor for cd=1
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<1,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , face_(org.face_)
      , isCopy_(org.isCopy_+1)
      , nrOfEdges_(org.nrOfEdges_)
      , item_(org.item_)
      , elem_(org.elem_)
      , iter_ ( org.iter_ )
  {

    indexList = new double[nrOfEdges_];
    for (int i = 0; i < nrOfEdges_; ++i)
      indexList[i] = org.indexList[i];

    // don't copy a copy of a copy of a copy of a copy
    assert( org.isCopy_ < 3 );
  }

  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLevelIterator<1, pitype, GridImp> :: increment () {

    if(endIter_)
      return ;
    IteratorType & iter = iter_.operator *();
    assert(face_>=0);

    int goNext = 1;
    item_ = &iter->getitem();
    while (face_ < 3) {
      int idx = item_->edge_idx(face_);
      if(!indexList[idx]) {
        indexList[idx]=1;
        goNext = 0;
        break;
      }
      ++face_;
    }

    if (goNext) {
      assert(face_==3);
      iter->next();
      if(iter->done()) {
        endIter_ = true;
        face_= 0;
        this->done();
        return ;
      }
      face_=0;
      item_ = &iter->getitem();
      this->updateEntityPointer(item_, face_, level_);
      increment();
      return;
    }

    if(iter->done()) {
      endIter_ = true;
      face_= 0;
      this->done();
      return ;
    }
    item_ = &iter->getitem();
    this->updateEntityPointer(item_, face_, level_);
    ++face_;
  }


  //********************************************************************
  //
  //  --ALU2dGridHierarchicIterator
  //  --HierarchicIterator
  //
  //********************************************************************


  //! the normal Constructor
  template<class GridImp>
  inline ALU2dGridHierarchicIterator<GridImp> ::
  ALU2dGridHierarchicIterator(const GridImp &grid, const HElementType & elem, int maxlevel, bool end) :
    ALU2dGridEntityPointer<0,GridImp> (grid)
    , elem_(&elem)
    , maxlevel_(maxlevel)
    , endIter_(end) {

    if (!end)
    {
      HElementType * item = const_cast<HElementType *> (elem_->down());
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


  //! the normal Constructor
  template<class GridImp>
  inline ALU2dGridHierarchicIterator<GridImp> ::
  ALU2dGridHierarchicIterator(const ALU2dGridHierarchicIterator<GridImp> &org) :
    ALU2dGridEntityPointer<0,GridImp> (org)
    , elem_ (org.elem_)
    , maxlevel_(org.maxlevel_)
    , endIter_(org.endIter_)  {}


  template <class GridImp>
  inline typename ALU2dGridHierarchicIterator<GridImp>::HElementType * ALU2dGridHierarchicIterator<GridImp>::
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

    nextelem = oldelem->father();
    if(nextelem == elem_) return 0;

    while( !nextelem->next() )
    {
      nextelem = nextelem->father();
      if(nextelem == elem_) return 0;
    }

    if(nextelem) nextelem = nextelem->next();

    return nextelem;
  }


  //! increment, go to next entity
  template<class GridImp>
  inline void ALU2dGridHierarchicIterator<GridImp> :: increment() {

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


} //end namespace Dune
#endif
