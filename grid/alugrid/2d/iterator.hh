// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDITERATOR_HH
#define DUNE_ALU2DGRIDITERATOR_HH

// System includes

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>

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
  class ALU2dGridIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator;
  template<int dim, int dimworld>
  class ALU2dGrid;


  //**********************************************************************
  //
  // --ALU2dGridIntersectionIterator
  // --IntersectionIterator
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, wh
     a neighbor is an entity of codimension 0 which has a common entity of codimens
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number o
     of an element!
   */
  template<class GridImp>
  class ALU2dGridIntersectionIterator
  //: public
  //  IntersectionIteratorDefaultImplementation<GridImp,
  //  ALU2dGridIntersectionIterator>
  {
    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    friend class IntersectionIteratorWrapper<GridImp>;
  public:
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

    //! neighbours
    //! The default Constructor , createing an empty ALU2dGridIntersectionIterator
    ALU2dGridIntersectionIterator(const GridImp & grid, int wLevel);

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU2dGridIntersectionIterator(const GridImp & grid,
                                  const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridIntersectionIterator(const ALU2dGridIntersectionIterator<GridImp> & org);

    //! The copy constructor
    void assign (const ALU2dGridIntersectionIterator<GridImp> & org);

    //! check whether entities are the same or whether iterator is done
    bool equals (const ALU2dGridIntersectionIterator<GridImp> & i) const;

    //! increment iterator
    void increment ();

    //! return level of inside(entity)
    int level () const;

    //! return true if intersection is with boundary
    bool boundary() const;

    int boundaryId () const;

    //! return true if intersection is with neighbor on this level
    bool neighbor () const;

    //! return true if intersection is with neighbor on this level
    bool leafNeighbor () const;

    //! return true if intersection is with neighbor on this level
    bool levelNeighbor () const;

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

    //! set interator to end iterator
    void done () ;

    //! reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    //! reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    // the local geometries
    mutable GeometryObject intersectionGlobal_;
    mutable GeometryObject intersectionSelfLocal_;
    mutable GeometryObject intersectionNeighborLocal_;

    // reference to grid
    const GridImp & grid_;

    //! current element from which we started the intersection iterator
    mutable HElementType* item_;
    mutable HElementType* neigh_;
    //mutable ALU2DSPACE Thinelement * nb_;

    mutable int nFaces_;
    mutable int walkLevel_;
    mutable int index_;

    mutable bool generatedGlobalGeometry_;
    mutable bool generatedLocalGeometries_;

    // unit outer normal
    mutable NormalType unitOuterNormal_;
    mutable NormalType outerNormal_;

    // true if end iterator
    bool done_;
    mutable bool calledOnLeaf_;
    mutable bool visitedLeaf_;

  }; // end ALU2dGridIntersectionIterator


  //********************************************************************
  //
  //  --TreeIterator
  //
  //********************************************************************

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class TreeIterator :
    public ALU2dGridEntityPointer<cdim,GridImp> {

    enum { dim = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim = cdim };

    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<cdim,GridImp> EntityPointerType;
    typedef ALU2dGridEntity<cdim,dim,GridImp> EntityImp;

  public:
    //! type of entity we iterate (interface)
    typedef typename GridImp::template Codim<cdim>::Entity Entity;
    typedef typename Dune::ALU2dImplTraits::template Codim<cdim>::InterfaceType ElementType;

    //! Constructor called by LeafIterator
    TreeIterator(const GridImp & grid, bool end);

    //! Constructor called by LevelIterator
    TreeIterator(const GridImp & grid, int level, bool end);

    //! copy Constructor
    TreeIterator(const TreeIterator<cdim, pitype, GridImp> & org);

    //! prefix increment
    void increment ();

  private:

    //! do not allow assigment
    TreeIterator<cdim, pitype, GridImp> & operator = (const TreeIterator<cdim, pitype, GridImp> & org)  {
      return *this;
    }

    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int face_;
    //! true if iterator is already a copy
    int isCopy_;

    ElementType * elem_;

    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    typedef typename ALU2dGridSpace:: AutoPointer< IteratorType >  IteratorPointerType ;
    IteratorPointerType iter_;

    const bool isLeafIterator_;
  }; // end TreeIterator



  //********************************************************************
  //
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************

  //! Leaf iterator
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator :
    public LeafIteratorDefaultImplementation<cdim, pitype, GridImp, ALU2dGridLeafIterator>,
    public TreeIterator<cdim,pitype,GridImp>
  {
    enum { dim = GridImp :: dimension };
    enum { dimworld  = GridImp::dimensionworld };
    enum { codim = cdim };

    friend class ALU2dGrid < dim , dimworld >;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<dim,dimworld,GridImp>;

    typedef ALU2dGridEntityPointer<cdim,GridImp> EntityPointerType;

  public:
    //! Constructor
    ALU2dGridLeafIterator(const GridImp & grid, bool end);
    //! copy Constructor
    ALU2dGridLeafIterator(const ALU2dGridLeafIterator<cdim, pitype, GridImp> & org);

  private:
    //! do not allow assigment
    ALU2dGridLeafIterator<cdim, pitype, GridImp> & operator = (const ALU2dGridLeafIterator<cdim, pitype, GridImp> & org)  {
      return *this;
    }
  }; // end ALU2dGridLeafIterator


  //**********************************************************************
  //
  // --ALU2dGridLevelIterator
  // --LevelIterator
  //**********************************************************************

  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */

  template<int cd, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLevelIterator :
    public LevelIteratorDefaultImplementation <cd, pitype, GridImp, ALU2dGridLevelIterator>,
    public TreeIterator<cd,pitype,GridImp>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class ALU2dGridEntity<dim,dimworld,GridImp>;
    friend class ALU2dGridEntity<1,dimworld,GridImp>;
    friend class ALU2dGridEntity<0,dimworld,GridImp>;
    friend class ALU2dGrid < dim , dimworld >;

    typedef ALU2dGridEntityPointer<cd,GridImp> EntityPointerType;

  public:
    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ALU2dGridLevelIterator<cd,pitype,GridImp> & org);

  private:
    //! do not allow assigment
    ALU2dGridLevelIterator<cd, pitype, GridImp> & operator = (const ALU2dGridLevelIterator<cd, pitype,GridImp> & org)  {
      return *this;
    }
  }; // end ALU2dGridLevelIterator


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

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ALU2dGridLevelIterator<2,pitype,GridImp> & org);

    //! prefix increment
    void increment ();

  private:
    //! do not allow assigment
    ALU2dGridLevelIterator<codim, pitype, GridImp> & operator = (const ALU2dGridLeafIterator<codim, pitype, GridImp> & org)  {
      return *this;
    }

    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int face_;
    //! true if iterator is already a copy
    int isCopy_;
    int nrOfVertices_;
    double* indexList;

    HElementType * item_;
    ALU2DSPACE Vertex * vertex_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits::template Codim<0>::InterfaceType ElementType;

    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    typedef typename ALU2dGridSpace:: AutoPointer< IteratorType >  IteratorPointerType ;
    IteratorPointerType iter_;

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

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    ALU2dGridLevelIterator(const GridImp & grid, int level, bool end);

    //! copy constructor
    ALU2dGridLevelIterator(const ALU2dGridLevelIterator<1,pitype,GridImp> & org);

    //! prefix increment
    void increment ();

  private:
    //! do not allow assigment
    ALU2dGridLevelIterator<codim, pitype, GridImp> & operator = (const ALU2dGridLeafIterator<codim, pitype, GridImp> & org)  {
      return *this;
    }

    //! true if iterator is end iterator
    bool endIter_;
    //! actual level
    int level_;
    //! information for edges
    int face_;
    //! true if iterator is already a copy
    int isCopy_;
    int nrOfEdges_;
    double* indexList;

    HElementType * item_;

    //! type of entity we iterate (interface)
    typedef typename Dune::ALU2dImplTraits::template Codim<1>::InterfaceType ElementType;
    ElementType * elem_;
    typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;
    typedef typename ALU2dGridSpace:: AutoPointer< IteratorType >  IteratorPointerType ;
    IteratorPointerType iter_;

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
