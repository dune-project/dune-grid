// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_HH
#define DUNE_ALU2DGRID_INTERSECTION_HH

// System includes
#include <stack>
#include <utility>

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>

#include <dune/grid/alugrid/2d/entity.hh>

namespace Dune
{

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridGeometry;
  template <class LocalGeometry, class LocalGeometryImp>
  class ALU2DIntersectionGeometryStorage;
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
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
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
  class ALU2dGridIntersectionBase
  {
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

  public:
    typedef ALU2dGridIntersectionBase< GridImp > ImplementationType;
    //! type of the intersection
    typedef Dune::Intersection< GridImp, Dune::ALU2dGridIntersectionBase > Intersection;

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dim,GridImp> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimworld> NormalType;

    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointerImp;

    typedef MakeableInterfaceObject< Geometry > GeometryObject;
    typedef MakeableInterfaceObject< LocalGeometry > LocalGeometryObject;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType;

    // type of local geometry storage
    typedef ALU2DIntersectionGeometryStorage<LocalGeometry, LocalGeometryImp>
    LocalGeometryStorageType ;

    typedef ALU2dGridIntersectionBase<GridImp> ThisType;
    friend class LevelIntersectionIteratorWrapper<GridImp>;
    friend class LeafIntersectionIteratorWrapper<GridImp>;

    friend class IntersectionIteratorWrapper<GridImp,ThisType>;
  protected:
    struct impl
    {
      impl() : item_(0) , neigh_(0) , index_(0) , opposite_(0), isBoundary_(false), useOutside_(false) { }
      impl(const impl & org) : item_(org.item_) , neigh_(org.neigh_) , index_(org.index_) , opposite_(org.opposite_),
                               isBoundary_(org.isBoundary_), useOutside_(org.useOutside_) { }

      impl & operator = (const impl & org)
      {
        item_ = org.item_;
        neigh_ = org.neigh_;
        index_ = org.index_;
        opposite_ = org.opposite_;
        isBoundary_ = org.isBoundary_;
        useOutside_ = org.useOutside_;
        return *this;
      }

      // current element from which we started the intersection iterator
      mutable HElementType* item_;
      mutable HElementType* neigh_;
      mutable int index_;
      mutable int opposite_;
      mutable bool isBoundary_;
      mutable bool useOutside_;
    } current;

  public:
    //! The default Constructor , creating an empty ALU2dGridIntersectionIterator
    ALU2dGridIntersectionBase(const GridImp & grid, int wLevel);

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU2dGridIntersectionBase(const GridImp & grid, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridIntersectionBase(const ThisType & org);

    virtual ~ALU2dGridIntersectionBase() {}

    //! The copy constructor
    void assign (const ThisType & org);

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }

    //! check whether entities are the same or whether iterator is done
    bool equals (const ThisType & i) const;

    //! return level of inside(entity)
    int level () const;

    //! return true if intersection is with boundary
    bool boundary() const;

    //! return boundary type
    int boundaryId () const;

    //! return the boundary segment index
    size_t boundarySegmentIndex() const;

    //! return true if intersection is with neighbor on this level
    bool neighbor () const;

    //! return EntityPointer to the Entity on the inside of this intersection.
    EntityPointer inside() const;

    //! return EntityPointer to the Entity on the outside of this intersection.
    EntityPointer outside() const;

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const;

    //! local index of codim 1 entity in neighbor where intersection is contained in
    int indexInOutside () const;

    int twistInInside () const;
    int twistInOutside () const;

    // deprecated methods
    int twistInSelf () const { return twistInInside(); }
    // deprecated methods
    int twistInNeighbor () const { return twistInOutside(); }

    NormalType & outerNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;
    NormalType & integrationOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;
    NormalType unitOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const;

    const LocalGeometry &geometryInInside () const;
    const LocalGeometry &geometryInOutside () const;
    const Geometry &geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const;

  protected:
    virtual bool conforming() const = 0;

    //! return true if intersection is with boundary
    void checkValid () ;

    // set interator to end iterator
    void done () ;

    // invalidate status of internal objects
    void unsetUp2Date() ;

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    // reset IntersectionIterator to first neighbour
    virtual void setFirstItem(const HElementType & elem, int wLevel);

    // the local geometries
    mutable GeometryObject intersectionGlobal_;
    mutable LocalGeometryObject intersectionSelfLocal_;
    mutable LocalGeometryObject intersectionNeighborLocal_;

    // reference to grid
    const GridImp & grid_;
    const LocalGeometryStorageType& localGeomStorage_;

    mutable int nFaces_;
    mutable int walkLevel_;

    // unit outer normal
    mutable NormalType outerNormal_;

    // true if end iterator
    bool done_;
  }; // end ALU2dGridIntersectionBase




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
    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType ;
    friend class LevelIntersectionIteratorWrapper<GridImp>;

    typedef ALU2dGridLevelIntersectionIterator<GridImp> ThisType;
    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

    typedef std::pair<HElementType *, std::pair<int, bool> > IntersectionInfo;

  public:
    typedef ALUMemoryProvider< ThisType > StorageType;

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

    typedef typename GridImp::template Codim<0>::Entity Entity;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dim,GridImp> LocalGeometryImp;
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

  protected:
    bool isConform() const {
      return this->neighbor() ?
             this->current.neigh_ == this->current.item_->neighbour(this->current.index_) :
             true;
    }
  public:
    //! level is conforming when non-conform grid used
    //! otherwise might not be conform
    bool conforming () const
    {
      return ( this->grid_.nonConform() ) ?
             true : isConform();
    }


  private:
    void doIncrement ();

    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    void addNeighboursToStack();

    int getOppositeInFather(int nrInChild, int nrOfChild) const;
    int getOppositeInChild(int nrInFather, int nrOfChild) const;

    mutable std::stack<IntersectionInfo> neighbourStack_;
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
    typedef ALU2dGridLeafIntersectionIterator<GridImp> ThisType;

    friend class LeafIntersectionIteratorWrapper<GridImp>;
    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

  public:
    typedef ALUMemoryProvider< ThisType > StorageType;

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

  private:
    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::ThinelementType ThinelementType;

    typedef std::pair< ThinelementType *, std::pair<int, bool> > IntersectionInfo;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dim,GridImp> EntityImp;
    typedef ALU2dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef ALU2dGridGeometry<dim-1,dim,GridImp> LocalGeometryImp;
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

    //! leaf is conforming, when conform grid version used
    bool conforming () const
    {
      return ( this->grid_.nonConform() ) ?
             ((this->neighbor()) ? (this->current.item_->level() == this->current.neigh_->level()) : true ) :
             true;
    }

  private:
    void doIncrement ();
    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    std::stack<IntersectionInfo> nbStack_;


  }; // end ALU2dGridLeafIntersectionIterator

} // end namespace Dune

#include "intersection_imp.cc"

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_HH
