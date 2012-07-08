// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_HH
#define DUNE_ALU2DGRID_INTERSECTION_HH

// System includes
#include <stack>
#include <utility>

// Dune includes
#include <dune/grid/common/grid.hh>

#include <dune/grid/alugrid/common/intersectioniteratorwrapper.hh>
#include <dune/grid/alugrid/2d/alu2dinclude.hh>
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



  // ALU2DIntersectionGeometryStorage
  // --------------------------------

  template< class LocalGeometryImpl >
  class ALU2DIntersectionGeometryStorage
  {
    typedef ALU2DIntersectionGeometryStorage< LocalGeometryImpl > ThisType;

    // one geometry for each face and twist 0 and 1
    LocalGeometryImpl geoms_[ 2 ][ 4 ][ 2 ];
    //std::vector< LocalGeometryImpl > geoms_[ 2 ][ 4 ];

  private:
    ALU2DIntersectionGeometryStorage ();

  public:
    // return reference to local geometry
    const LocalGeometryImpl &localGeom ( const int aluFace, const int twist, const int corners ) const
    {
      assert( corners == 3 || corners == 4 );
      assert( 0 <= aluFace && aluFace < corners );
      assert( twist == 0 || twist == 1 );
      return geoms_[ corners-3 ][ aluFace ][ twist ];
    }

    // return static instance
    static const ThisType &instance ()
    {
      static const ThisType geomStorage;
      return geomStorage;
    }
  };



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

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;

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
    typedef FieldVector< alu2d_ctype, dimworld > NormalType;
    typedef FieldVector< alu2d_ctype, dim-1 > LocalCoordinate;

    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointerImp;

    typedef typename ALU2dImplTraits< dimworld, eltype >::ThinelementType ThinelementType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::HBndElType HBndElType;

    // type of local geometry storage
    typedef ALU2DIntersectionGeometryStorage< LocalGeometryImpl > LocalGeometryStorageType;

    typedef ALU2dGridIntersectionBase<GridImp> ThisType;
    friend class LevelIntersectionIteratorWrapper<GridImp>;
    friend class LeafIntersectionIteratorWrapper<GridImp>;

    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

  protected:
    struct impl
    {
      explicit impl( HElementType *inside = 0 )
        : index_(0), useOutside_(false)
      {
        setInside( inside );
        setOutside( 0, -1 );
      }

      HElementType *inside () const
      {
        return inside_;
      }

      int nFaces () const
      {
        return nFaces_;
      }

      bool isBoundary () const
      {
        assert( inside() && (index_ < nFaces()) && inside()->neighbour( index_ ) );
        return inside()->neighbour( index_ )->thinis( ThinelementType::bndel_like );
      }

      HBndElType *boundary () const
      {
        assert( isBoundary() );
        return (HBndElType *)inside()->neighbour( index_ );
      }

      HElementType *outside () const
      {
        return outside_;
      }

      int opposite () const
      {
        return opposite_;
      }

      void setInside ( HElementType *inside )
      {
        inside_ = inside;
        nFaces_ = (inside != 0 ? inside->numfaces() : 0);
      }

      void setOutside ( HElementType *outside, int opposite )
      {
        outside_ = outside;
        opposite_ = opposite;
      }

    private:
      // current element from which we started the intersection iterator
      HElementType *inside_;
      HElementType *outside_;
      int nFaces_;
      int opposite_;

    public:
      mutable int index_;
      mutable bool useOutside_;
    } current;

  public:
    //! The default Constructor , creating an empty ALU2dGridIntersectionIterator
    ALU2dGridIntersectionBase(const FactoryType& factory, int wLevel);

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

    //! \deprecated Use twistInInside() instead. This method will be removed after Dune 2.3.
    int twistInSelf () const
    DUNE_DEPRECATED_MSG("Use twistInInside() instead.")
    { return twistInInside(); }

    //! \deprecated Use twistInOutside() instead. This method will be removed after Dune 2.3.
    int twistInNeighbor () const
    DUNE_DEPRECATED_MSG("Use twistInOutside() instead.")
    { return twistInOutside(); }

    NormalType outerNormal ( const LocalCoordinate &local ) const;
    NormalType integrationOuterNormal ( const LocalCoordinate &local ) const;
    NormalType unitOuterNormal ( const LocalCoordinate &local ) const;

    LocalGeometry geometryInInside () const;
    LocalGeometry geometryInOutside () const;
    Geometry geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const;

  protected:
    const GridImp& grid() const { return factory_.grid(); }

    virtual bool conforming() const = 0;

    //! return true if intersection is with boundary
    void checkValid () ;

    // set interator to end iterator
    void done ( const HElementType *inside );
    void done ( const EntityImp &en ) { done( &en.getItem() ); }

    // invalidate status of internal objects
    void unsetUp2Date() ;

    // reset IntersectionIterator to first neighbour
    void first ( const EntityImp &en, int wLevel );

    // reset IntersectionIterator to first neighbour
    virtual void setFirstItem ( const HElementType &elem, int wLevel ) = 0;

    // the local geometries
    mutable GeometryImpl intersectionGlobal_;
    mutable LocalGeometryImpl intersectionSelfLocal_;
    mutable LocalGeometryImpl intersectionNeighborLocal_;

    // reference to factory
    const FactoryType& factory_;
    const LocalGeometryStorageType &localGeomStorage_;

    mutable int walkLevel_;
  }; // end ALU2dGridIntersectionBase




  //**********************************************************************
  //
  // --ALU2dGridLevelIntersectionIterator
  // --IntersectionLevelIterator
  //**********************************************************************

  template< class GridImp >
  class ALU2dGridLevelIntersectionIterator
    : public ALU2dGridIntersectionBase< GridImp >
  {
    typedef ALU2dGridLevelIntersectionIterator< GridImp > ThisType;
    typedef ALU2dGridIntersectionBase< GridImp > BaseType;

    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    typedef typename ALU2dImplTraits< dimworld, eltype >::ThinelementType ThinelementType;
    typedef typename BaseType::HElementType HElementType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::HBndElType HBndElType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::PeriodicBndElType PeriodicBndElType;

    friend class LevelIntersectionIteratorWrapper<GridImp>;

    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

    typedef std::pair< HElementType *, int > IntersectionInfo;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;
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
    ALU2dGridLevelIntersectionIterator(const FactoryType& factory, int wLevel);

    //! The default Constructor , level tells on which level we want neighbours
    ALU2dGridLevelIntersectionIterator(const FactoryType& factory, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridLevelIntersectionIterator(const ALU2dGridLevelIntersectionIterator<GridImp> & org);

    void assign (const ALU2dGridLevelIntersectionIterator<GridImp> & org);

    //! increment iterator
    void increment ();

  public:
    //! level is conforming when non-conform grid used
    //! otherwise might not be conform
    bool conforming () const
    {
      return (this->grid().nonConform() || isConform());
    }

  protected:
    bool isConform() const
    {
      return (!current.outside() || (current.outside() == current.inside()->neighbour( current.index_ )));
    }

  private:
    void doIncrement ();

    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    void addNeighboursToStack();

    static int getOppositeInFather ( int nrInChild, int nrOfChild );
    static int getOppositeInChild ( int nrInFather, int nrOfChild );

    void setupIntersection ();

  protected:
    using BaseType::done;

    using BaseType::current;
    using BaseType::walkLevel_;

  private:
    mutable std::stack<IntersectionInfo> nbStack_;
  }; // end ALU2dGridLevelIntersectionIterator



  //********************************************************************
  //
  //  --ALU2dGridLeafIntersectionIterator
  //
  //
  //********************************************************************

  template< class GridImp >
  class ALU2dGridLeafIntersectionIterator
    : public ALU2dGridIntersectionBase< GridImp >
  {
    typedef ALU2dGridLeafIntersectionIterator< GridImp > ThisType;
    typedef ALU2dGridIntersectionBase< GridImp > BaseType;

    friend class LeafIntersectionIteratorWrapper<GridImp>;
    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

    static const int dim = GridImp::dimension;
    static const int dimworld  = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

  public:
    typedef typename GridImp :: GridObjectFactoryType FactoryType;
    typedef ALUMemoryProvider< ThisType > StorageType;

    enum { dimension       = GridImp::dimension };
    enum { dimensionworld  = GridImp::dimensionworld };

  private:
    typedef typename ALU2dImplTraits< dimworld, eltype >::ThinelementType ThinelementType;
    typedef typename BaseType::HElementType HElementType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::HBndElType HBndElType;
    typedef typename ALU2dImplTraits< dimworld, eltype >::PeriodicBndElType PeriodicBndElType;

    typedef std::pair< HElementType *, int > IntersectionInfo;

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
    ALU2dGridLeafIntersectionIterator(const FactoryType& factory, int wLevel);

    //! The default Constructor , level tells on which level we want neighbours
    ALU2dGridLeafIntersectionIterator(const FactoryType& factory, const HElementType* el, int wLevel, bool end=true);

    //! The copy constructor
    ALU2dGridLeafIntersectionIterator(const ALU2dGridLeafIntersectionIterator<GridImp> & org);

    void assign (const ALU2dGridLeafIntersectionIterator<GridImp> & org);

    //! increment iterator
    void increment ();

  public:
    //! leaf is conforming, when conform grid version used
    bool conforming () const
    {
      return (!this->grid().nonConform() || isConform());
    }

  protected:
    bool isConform() const
    {
      return (!current.outside() || (current.outside()->level() == current.inside()->level()) );
    }

  private:
    void doIncrement ();
    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    void setupIntersection ();

  protected:
    using BaseType::done;

    using BaseType::current;
    using BaseType::walkLevel_;

  private:
    std::stack<IntersectionInfo> nbStack_;

  }; // end ALU2dGridLeafIntersectionIterator

} // end namespace Dune

#include "intersection_imp.cc"
#if COMPILE_ALU2DGRID_INLINE
  #include "intersection.cc"
#endif

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_HH
