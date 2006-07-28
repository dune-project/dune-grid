// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDITERATOR_HH
#define DUNE_ALU2DGRIDITERATOR_HH

// System includes

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>
#include <stack>
#include <utility>


// Local includes
#include "entity.hh"
#include "myautoptr.hh"

namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridGeometry;
  template<class GridImp>
  class ALU2dGridHierarchicIterator;
  template<class GridImp>
  class ALU2dGridIntersectionBase;
  template<class GridImp>
  class ALU2dGridLeafIntersectionIterator;
  template<class GridImp>
  class ALU2dGridLevelIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator;
  template<int dim, int dimworld>
  class ALU2dGrid;


  //**********************************************************************
  //
  // --ALU2dGridIntersectionBase
  // --IntersectionBase
  //**********************************************************************
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, wh
     a neighbor is an entity of codimension 0 which has a common entity of codimens
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number o
     of an element!
   */

  template<class GridImp>
  class ALU2dGridIntersectionBase : public IntersectionIteratorDefaultImplementation <GridImp, ALU2dGridLevelIntersectionIterator >
  {
  public:
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimworld> NormalType;

    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointerImp;

    typedef MakeableInterfaceObject< Geometry > GeometryObject;

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    friend class ALU2dGridLevelIntersectionIterator<GridImp>;
    friend class ALU2dGridLeafIntersectionIterator<GridImp>;

  protected:
    struct impl
    {
      impl() : item_(0) , neigh_(0) , index_(0) , opposite_(0) {}
      impl(const impl & org) : item_(org.item_) , neigh_(org.neigh_) , index_(org.index_) , opposite_(org.opposite_) {}

      impl & operator = (const impl & org)
      {
        item_ = org.item_;
        neigh_ = org.neigh_;
        index_ = org.index_;
        opposite_ = org.opposite_;
        return *this;
      }
      // current element from which we started the intersection iterator
      mutable HElementType* item_;
      mutable HElementType* neigh_;
      mutable int index_;
      mutable int opposite_;
    } current;

  public:
    //! The default Constructor , creating an empty ALU2dGridIntersectionIterator
    ALU2dGridIntersectionBase(const GridImp & grid, int wLevel);

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU2dGridIntersectionBase(const GridImp & grid, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridIntersectionBase(const ALU2dGridIntersectionBase<GridImp> & org);

    virtual ~ALU2dGridIntersectionBase() {}

    //! The copy constructor
    void assign (const ALU2dGridIntersectionBase<GridImp> & org);

    //! check whether entities are the same or whether iterator is done
    bool equals (const ALU2dGridIntersectionBase<GridImp> & i) const;

    //! increment iterator
    virtual void increment() { };

    //! return level of inside(entity)
    int level () const;

    //! return true if intersection is with boundary
    bool boundary() const;

    int boundaryId () const;

    //! return true if intersection is with neighbor on this level
    bool neighbor () const;

    //! return EntityPointer to the Entity on the inside of this intersection.
    EntityPointer inside() const;

    //! return EntityPointer to the Entity on the outside of this intersection.
    EntityPointer outside() const;

    //! local number of codim 1 entity in self where intersection is contained in
    int numberInSelf () const;

    //! local number of codim 1 entity in neighbor where intersection is contained in
    int numberInNeighbor () const;

    NormalType & outerNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;
    NormalType & integrationOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;
    NormalType & unitOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;

    const LocalGeometry & intersectionSelfLocal () const;
    const LocalGeometry & intersectionNeighborLocal () const;
    const Geometry & intersectionGlobal () const;

  private:

    // set interator to end iterator
    void done () ;
    // reset IntersectionIterator to first neighbour
    virtual void setFirstItem(const HElementType & elem, int wLevel);
    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    // the local geometries
    mutable GeometryObject intersectionGlobal_;
    mutable GeometryObject intersectionSelfLocal_;
    mutable GeometryObject intersectionNeighborLocal_;

    // reference to grid
    const GridImp & grid_;
    mutable int nFaces_;
    mutable int walkLevel_;
    mutable bool generatedGlobalGeometry_;
    mutable bool generatedLocalGeometries_;

    // unit outer normal
    mutable NormalType unitOuterNormal_;
    mutable NormalType outerNormal_;

    // true if end iterator
    bool done_;
  }; // end ALU2dGridIntersectionIterator




  //**********************************************************************
  //
  // --ALU2dGridLevelIntersectionIterator
  // --IntersectionLevelIterator
  //**********************************************************************

  template<class GridImp>
  class ALU2dGridLevelIntersectionIterator
    : public ALU2dGridIntersectionBase<GridImp>
      //,  public IntersectionIteratorDefaultImplementation <GridImp, ALU2dGridLevelIntersectionIterator >
  {

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    friend class LevelIntersectionIteratorWrapper<GridImp>;

    typedef ALU2dGridLevelIntersectionIterator<GridImp> ThisType;
  public:
    typedef ALUMemoryProvider< ThisType > StorageType;

    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

    typedef typename GridImp::template Codim<0>::Entity Entity;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimworld> NormalType;
    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointer;

    typedef MakeableInterfaceObject< Geometry > GeometryObject;

    //! The default Constructor , creating an empty ALU2dGridIntersectionIterator
    ALU2dGridLevelIntersectionIterator(const GridImp & grid, int wLevel);

    //! The default Constructor , level tells on which level we want neighbours
    ALU2dGridLevelIntersectionIterator(const GridImp & grid, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridLevelIntersectionIterator(const ALU2dGridLevelIntersectionIterator<GridImp> & org);

    void assign (const ALU2dGridLevelIntersectionIterator<GridImp> & org);

    //! increment iterator
    void increment ();

  private:
    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    void addNeighboursToStack();

    int getOppositeInFather(int nrInChild, int nrOfChild) const;
    int getOppositeInChild(int nrInFather, int nrOfChild) const;

    mutable std::stack<std::pair<HElementType *, int> > neighbourStack_;
  }; // end ALU2dGridLevelIntersectionIterator


  //********************************************************************
  //
  //  --ALU2dGridLeafIntersectionIterator
  //
  //
  //********************************************************************

  template<class GridImp>
  class ALU2dGridLeafIntersectionIterator
    : public ALU2dGridIntersectionBase<GridImp>
      //, public IntersectionIteratorDefaultImplementation<GridImp, ALU2dGridLevelIntersectionIterator>
  {
    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    friend class LeafIntersectionIteratorWrapper<GridImp>;

    typedef ALU2dGridLeafIntersectionIterator<GridImp> ThisType;
  public:
    typedef ALUMemoryProvider< ThisType > StorageType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimworld> NormalType;
    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointer;

    typedef MakeableInterfaceObject< Geometry > GeometryObject;

    //! The default Constructor , createing an empty ALU2dGridIntersectionIterator
    ALU2dGridLeafIntersectionIterator(const GridImp & grid, int wLevel);

    //! The default Constructor , level tells on which level we want neighbours
    ALU2dGridLeafIntersectionIterator(const GridImp & grid, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridLeafIntersectionIterator(const ALU2dGridLeafIntersectionIterator<GridImp> & org);

    void assign (const ALU2dGridLeafIntersectionIterator<GridImp> & org);

    //! increment iterator
    void increment ();

  private:
    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    std::stack<ALU2DSPACE Thinelement*> nbStack_;


  }; // end ALU2dGridLeafIntersectionIterator



  //********************************************************************
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator :
    public ALU2dGridEntityPointer<cdim,GridImp>,
    public LeafIteratorDefaultImplementation<cdim, pitype, GridImp, ALU2dGridLeafIterator>
  {
    enum { dim = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim = cdim };

    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<cdim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<cdim,dim,GridImp> EntityImp;

    typedef ALU2dGridLeafIterator<cdim, pitype, GridImp> ThisType;

  public:
    //! type of entity we iterate (interface)
    typedef typename GridImp::template Codim<cdim>::Entity Entity;
    typedef typename Dune::ALU2dImplTraits::template Codim<cdim>::InterfaceType ElementType;

    //! Constructor called by LeafIterator
    ALU2dGridLeafIterator(const GridImp & grid, bool end);

    //! copy Constructor
    ALU2dGridLeafIterator(const ThisType & org);

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int face_;

    ElementType * elem_;

    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

    // Listwalkptr, behaves like a proxy for Leafwalk and Levelwalk Ptrs
    IteratorType iter_;

  }; // end ALU2dGridLeafIterator



  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=0
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<0, pitype, GridImp> :
    public LevelIteratorDefaultImplementation <0, pitype, GridImp, ALU2dGridLevelIterator>,
    public ALU2dGridEntityPointer<0,GridImp>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim  = 0 };

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    typedef ALU2dGridLevelIterator<0,pitype,GridImp> ThisType;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    //! prefix increment
    void increment ();

    //! do not allow assigment
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;

    HElementType * item_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits::template Codim<0>::InterfaceType ElementType;
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    IteratorType iter_;

  };

  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=1
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<1, pitype, GridImp> :
    public LevelIteratorDefaultImplementation <1, pitype, GridImp, ALU2dGridLevelIterator>,
    public ALU2dGridEntityPointer<1,GridImp>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim  = 1 };

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    typedef ALU2dGridLevelIterator<1,pitype,GridImp> ThisType;
  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    ~ALU2dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int myFace_;

    // current item
    HElementType * item_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits::template Codim<1>::InterfaceType ElementType;
    ElementType * elem_;
    // Listwalkptr is a proxy for iterator pointers
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

    IteratorType iter_;

    ALU2dGridMarkerVector & marker_;
  };

  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator, specialized for cd=2
  //**********************************************************************

  template<PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator<2, pitype, GridImp> :
    public LevelIteratorDefaultImplementation <2, pitype, GridImp, ALU2dGridLevelIterator>,
    public ALU2dGridEntityPointer<2,GridImp>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim  = 2 };

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<codim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<codim,dim,GridImp> EntityImp;

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;
    typedef ALU2dGridLevelIterator<2,pitype,GridImp> ThisType;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ThisType & org);

    ~ALU2dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! assigment of iterator
    ThisType & operator = (const ThisType & org);

  private:
    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int myFace_;
    //! true if iterator is already a copy

    HElementType * item_;
    ALU2DSPACE Vertex * vertex_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits::template Codim<0>::InterfaceType ElementType;

    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    IteratorType iter_;
    ALU2dGridMarkerVector & marker_;
  };

  //***************************************************************
  //
  // - HierarchicIteraror
  // --HierarchicIterator
  //***************************************************************

  //! Hierarichic Iterator of ALU2dGrid
  template<class GridImp>
  class ALU2dGridHierarchicIterator :
    public ALU2dGridEntityPointer<0,GridImp> ,
    public HierarchicIteratorDefaultImplementation <GridImp,ALU2dGridHierarchicIterator>
  {
    enum { dim = GridImp::dimension };
    // my type
    typedef ALU2dGridHierarchicIterator<GridImp> ThisType;

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

  public:
    //! type of entities we iterate
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! type of coordinates, i.e. double
    typedef typename GridImp::ctype ctype;
    //! tpye of entity implementation
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;

    //! the normal Constructor
    ALU2dGridHierarchicIterator(const GridImp &grid, const HElementType & elem, int maxlevel, bool end=false);

    //! the normal Constructor
    ALU2dGridHierarchicIterator(const ALU2dGridHierarchicIterator<GridImp> &org);

    //! increment, go to next entity
    void increment();

    //! the assignment operator
    ThisType & operator = (const ALU2dGridHierarchicIterator<GridImp> &org)
    {
      ALU2dGridEntityPointer<0,GridImp> :: operator = (org);
      elem_ = org.elem_;
      maxlevel_= org.maxlevel_;
      endIter_ = org.endIter_;
      return *this;
    };

  private:

    //! go to next valid element
    HElementType * goNextElement (HElementType * oldEl);

    //! element from where we started
    const HElementType * elem_;

    //! maximal level to go down
    int maxlevel_;
    //! true if iterator is end iterator
    bool endIter_;
  }; // end ALU2dHierarchicIterator

} // end namespace Dune

#include "iterator_imp.cc"

#endif
