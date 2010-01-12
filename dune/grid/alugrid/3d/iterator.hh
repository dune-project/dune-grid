// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDITERATOR_HH
#define DUNE_ALU3DGRIDITERATOR_HH

// System includes

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"
#include "faceutility.hh"
#include "alu3diterators.hh"
#include "memory.hh"

namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU3dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU3dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU3dGridGeometry;
  template<class GridImp>
  class ALU3dGridHierarchicIterator;
  template<class GridImp>
  class ALU3dGridIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLeafIterator;
  template<int dim, int dimworld, ALU3dGridElementType elType>
  class ALU3dGrid;
  template <ALU3dGridElementType type>
  class ALU3dGridFaceInfo;
  template <ALU3dGridElementType elType>
  class ALU3dGridGeometricFaceInfo;

  //**********************************************************************
  //
  // --ALU3dGridIntersectionIterator
  // --IntersectionIterator
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, wh
     a neighbor is an entity of codimension 0 which has a common entity of codimens
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number o
     of an element!
   */
  template<class GridImp>
  class ALU3dGridIntersectionIterator
  //: public IntersectionIteratorDefaultImplementation <GridImp,ALU3dGridIntersectionIterator>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef ALU3dImplTraits<GridImp::elementType> ImplTraits;

    typedef typename ImplTraits::HElementType HElementType ;
    typedef typename ImplTraits::HBndSegType HBndSegType;
    typedef typename ImplTraits::GEOElementType GEOElementType;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::GEOFaceType GEOFaceType;
    typedef typename ImplTraits::NeighbourPairType NeighbourPairType;
    typedef typename ImplTraits::PLLBndFaceType PLLBndFaceType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

    typedef ALU3DSPACE GEOTetraElementType GEOTetraElementType;
    typedef ALU3DSPACE GEOHexaElementType GEOHexaElementType;
    typedef ALU3DSPACE BNDFace3Type GEOTriangleBndType;
    typedef ALU3DSPACE BNDFace4Type GEOQuadBndType;

    typedef ALU3dGridFaceInfo<GridImp::elementType> FaceInfoType;
    typedef typename std::auto_ptr<FaceInfoType> FaceInfoPointer;

    typedef typename SelectType<
        is_same<Int2Type<tetra>, Int2Type<GridImp::elementType> >::value,
        ALU3dGridGeometricFaceInfoTetra,
        ALU3dGridGeometricFaceInfoHexa
        >::Type GeometryInfoType;

    typedef ElementTopologyMapping<GridImp::elementType> ElementTopo;
    typedef FaceTopologyMapping<GridImp::elementType> FaceTopo;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    enum { numVerticesPerFace =
             EntityCount<GridImp::elementType>::numVerticesPerFace };
    enum { numVertices = EntityCount<GridImp::elementType>::numVertices };

    typedef ALU3dGridIntersectionIterator<GridImp> ThisType;

    friend class ALU3dGridEntity<0,dim,GridImp>;
    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

  protected:
    enum IntersectionIteratorType { IntersectionLeaf , IntersectionLevel, IntersectionBoth };

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef ALU3dGridIntersectionIterator< GridImp > ImplementationType;
    //! type of the intersection
    typedef Dune::Intersection< GridImp, Dune::ALU3dGridIntersectionIterator > Intersection;

    typedef ALU3dGridGeometry<dim-1,dimworld,GridImp> GeometryImp;
    typedef MakeableInterfaceObject<Geometry> GeometryObject;

    typedef FieldVector<alu3d_ctype, dimworld> NormalType;
    typedef ALU3dGridEntityPointer<0,GridImp> EntityPointer;

    typedef ALUMemoryProvider< ThisType > StorageType;

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU3dGridIntersectionIterator(const GridImp & grid,
                                  HElementType *el,
                                  int wLevel,bool end=false);

    ALU3dGridIntersectionIterator(const GridImp & grid,int wLevel);

    //! The copy constructor
    ALU3dGridIntersectionIterator(const ALU3dGridIntersectionIterator<GridImp> & org);

    //! assignment of iterators
    void assign(const ALU3dGridIntersectionIterator<GridImp> & org);

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }

    //! The copy constructor
    bool equals (const ALU3dGridIntersectionIterator<GridImp> & i) const;

    //! increment iterator
    void increment ();

    //! access neighbor
    EntityPointer outside() const;

    //! access entity where iteration started
    EntityPointer inside() const;

    //! return true if intersection is with boundary.
    bool boundary () const;

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const;

    //! return true if across the edge an neighbor on this level exists
    bool levelNeighbor () const;

    //! return true if across the edge an neighbor on leaf level exists
    bool leafNeighbor () const;

    //! return information about the Boundary
    int boundaryId () const;

    //! return the boundary segment index
    size_t boundarySegmentIndex() const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const;

    //! intersection of codimension 1 of this neighbor with element where
    //!  iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where
    //! iteration started.
    const Geometry &geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const;

    //! local index of codim 1 entity in self where intersection is contained
    //!  in
    int indexInInside () const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const;

    //! local index of codim 1 entity in neighbor where intersection is
    //! contained
    int indexInOutside () const;

    //! returns twist of face compared to inner element
    int twistInSelf() const { return twistInInside(); }

    //! returns twist of face compared to outer element
    int twistInNeighbor() const { return twistInOutside(); }

    //! returns twist of face compared to inner element
    int twistInInside() const;

    //! returns twist of face compared to outer element
    int twistInOutside() const;

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & unitOuterNormal (const FieldVector<alu3d_ctype, dim-1>& local) const ;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & outerNormal (const FieldVector<alu3d_ctype, dim-1>& local) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    NormalType & integrationOuterNormal (const FieldVector<alu3d_ctype, dim-1>& local) const;

    //! return level of iterator
    int level () const;

    //! return true if intersection is conforming
    bool conforming () const
    {
      return (connector_.conformanceState() == FaceInfoType::CONFORMING);
    }

  protected:
    // set interator to end iterator
    void done () ;

    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    void setInteriorItem(const HElementType  & elem,
                         const PLLBndFaceType& bnd, int wLevel);

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    // set new face
    void setNewFace(const GEOFaceType& newFace);

  private:
    // set new face (only LeafIntersectionIterator)
    void setGhostFace(const GEOFaceType& newFace);

  protected:
    // generate local geometries
    void buildLocalGeometries() const;

    // get the face corresponding to the index
    const typename ALU3dImplTraits<tetra>::GEOFaceType*
    getFace(const GEOTriangleBndType & bnd, int index) const;

    // get the face corresponding to the index
    const typename ALU3dImplTraits<hexa>::GEOFaceType*
    getFace(const GEOQuadBndType & bnd, int index) const;

    // get the face corresponding to the index
    const typename ALU3dImplTraits<tetra>::GEOFaceType*
    getFace(const GEOTetraElementType & elem, int index) const;

    const typename ALU3dImplTraits<hexa>::GEOFaceType*
    getFace(const GEOHexaElementType & elem, int index) const;

    //! structure containing the topological and geometrical information about
    //! the face which the iterator points to
    mutable FaceInfoType connector_;
    mutable GeometryInfoType geoProvider_; // need to initialise

    // reference to grid
    const GridImp & grid_;

    //! current element from which we started the intersection iterator
    const IMPLElementType* item_;

    //! current pointer to ghost face if iterator was started from ghost element
    const PLLBndFaceType* ghost_;

    mutable int innerLevel_;
    mutable int index_;

    mutable GeometryObject intersectionGlobal_;
    mutable GeometryImp &  intersectionGlobalImp_;
    mutable GeometryObject intersectionSelfLocal_;
    mutable GeometryImp &  intersectionSelfLocalImp_;
    mutable GeometryObject intersectionNeighborLocal_;
    mutable GeometryImp &  intersectionNeighborLocalImp_;

    // unit outer normal
    mutable NormalType unitOuterNormal_;

    // true if end iterator
    bool done_;
  };

  template<class GridImp>
  class ALU3dGridLevelIntersectionIterator :
    public ALU3dGridIntersectionIterator<GridImp>
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef ALU3dImplTraits<GridImp::elementType> ImplTraits;
    typedef ALU3DSPACE HElementType HElementType ;
    typedef typename ImplTraits::GEOElementType GEOElementType;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::GEOFaceType GEOFaceType;
    typedef typename ImplTraits::NeighbourPairType NeighbourPairType;
    typedef typename ImplTraits::PLLBndFaceType PLLBndFaceType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

    typedef ALU3dGridFaceInfo<GridImp::elementType> FaceInfoType;
    typedef typename std::auto_ptr<FaceInfoType> FaceInfoPointer;

    typedef typename SelectType<
        is_same<Int2Type<tetra>, Int2Type<GridImp::elementType> >::value,
        ALU3dGridGeometricFaceInfoTetra,
        ALU3dGridGeometricFaceInfoHexa
        >::Type GeometryInfoType;

    typedef ElementTopologyMapping<GridImp::elementType> ElementTopo;
    typedef FaceTopologyMapping<GridImp::elementType> FaceTopo;

    enum { numFaces = EntityCount<GridImp::elementType>::numFaces };
    enum { numVerticesPerFace =
             EntityCount<GridImp::elementType>::numVerticesPerFace };
    enum { numVertices = EntityCount<GridImp::elementType>::numVertices };

    typedef ALU3dGridLevelIntersectionIterator<GridImp> ThisType;

    friend class ALU3dGridEntity<0,dim,GridImp>;
    friend class IntersectionIteratorWrapper<GridImp,ThisType>;

  public:
    //typedef Dune :: Intersection< const GridImp, ThisType >
    typedef ALUMemoryProvider< ThisType > StorageType;

    //! The default Constructor , level tells on which level we want
    //! neighbours
    ALU3dGridLevelIntersectionIterator(const GridImp & grid,
                                       HElementType *el,
                                       int wLevel,bool end=false);

    ALU3dGridLevelIntersectionIterator(const GridImp & grid,int wLevel);

    //! The copy constructor
    ALU3dGridLevelIntersectionIterator(const ThisType & org);

    //! assignment of iterators
    void assign(const ThisType & org);

    //! increment iterator
    void increment ();

    // reset IntersectionIterator to first neighbour
    template <class EntityType>
    void first(const EntityType & en, int wLevel);

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const;

    //! return true if across the edge an neighbor on this level exists
    bool levelNeighbor () const;

    //! return true if across the edge an neighbor on leaf level exists
    bool leafNeighbor () const;

    //! return true if intersection is conforming
    bool conforming () const
    {
      assert( !neighbor() || this->connector_.conformanceState() == FaceInfoType::CONFORMING );
      return true;
    }
  private:
    // set new face
    void setNewFace(const GEOFaceType& newFace);

    // reset IntersectionIterator to first neighbour
    void setFirstItem(const HElementType & elem, int wLevel);

    // reset IntersectionIterator to first neighbour
    void setInteriorItem(const HElementType  & elem,
                         const PLLBndFaceType& bnd, int wLevel);

    bool levelNeighbor_;
    bool isLeafItem_;
  };

  //////////////////////////////////////////////////////////////////////////////
  //
  //  --IterationImpl
  //
  //////////////////////////////////////////////////////////////////////////////
  template <class InternalIteratorType >
  struct ALU3dGridTreeIterator
  {
    typedef typename InternalIteratorType :: val_t val_t;

    // here the items level will do
    template <class GridImp, int codim>
    struct GetLevel
    {
      template <class ItemType>
      static int getLevel(const GridImp & grid, const ItemType & item, int level )
      {
        assert( & item );
        return (level < 0) ? item.level() : level;
      }
    };

    // level is not needed for codim = 0
    template <class GridImp>
    struct GetLevel<GridImp,0>
    {
      template <class ItemType>
      static int getLevel(const GridImp & grid, const ItemType & item, int level )
      {
        return level;
      }
    };

    template <class GridImp>
    struct GetLevel<GridImp,3>
    {
      template <class ItemType>
      static int getLevel(const GridImp & grid, const ItemType & item, int level)
      {
        return (level < 0) ? grid.getLevelOfLeafVertex(item) : level;
      }
    };

  protected:
    // set iterator to first item
    template <class GridImp, class IteratorImp>
    void firstItem(const GridImp & grid, IteratorImp & it, int level )
    {
      InternalIteratorType & iter = it.internalIterator();
      iter.first();
      if( ! iter.done() )
      {
        assert( iter.size() > 0 );
        setItem(grid,it,iter,level);
      }
      else
      {
        it.removeIter();
      }
    }

    // set the iterators entity to actual item
    template <class GridImp, class IteratorImp>
    void setItem (const GridImp & grid, IteratorImp & it, InternalIteratorType & iter, int level)
    {
      enum { codim = IteratorImp :: codimension };
      val_t & item = iter.item();
      assert( item.first || item.second );
      if( item.first )
      {
        it.updateEntityPointer( item.first ,
                                GetLevel<GridImp,codim>::getLevel(grid, *(item.first) , level) );
      }
      else
        it.updateGhostPointer( *item.second );
    }

    // increment iterator
    template <class GridImp, class IteratorImp>
    void incrementIterator(const GridImp & grid, IteratorImp & it, int level)
    {
      // if iter_ is zero, then end iterator
      InternalIteratorType & iter = it.internalIterator();

      iter.next();

      if(iter.done())
      {
        it.removeIter();
        return ;
      }

      setItem(grid,it,iter,level);
      return ;
    }
  };

  //**********************************************************************
  //
  // --ALU3dGridLevelIterator
  // --LevelIterator
  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template<int cd, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLevelIterator
    : public ALU3dGridEntityPointer <cd,GridImp>,
      // public LevelIteratorDefaultImplementation <cd,pitype,GridImp,ALU3dGridLevelIterator> ,
      public ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> >
  {
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class ALU3dGridEntity<3,dim,GridImp>;
    friend class ALU3dGridEntity<2,dim,GridImp>;
    friend class ALU3dGridEntity<1,dim,GridImp>;
    friend class ALU3dGridEntity<0,dim,GridImp>;
    friend class ALU3dGrid < dim , dimworld, GridImp::elementType >;

    friend class ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> > ;
  public:
    typedef typename GridImp::template Codim<cd>::Entity Entity;
    typedef ALU3DSPACE VertexListType VertexListType;

    //! typedef of my type
    typedef ALU3dGridLevelIterator<cd,pitype,GridImp> ThisType;
    // the wrapper for the original iterator of the ALU3dGrid
    typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper<cd,pitype> IteratorType;
    typedef IteratorType InternalIteratorType;
    typedef typename ALU3DSPACE IteratorElType<cd>::val_t val_t;

    //! Constructor for begin iterator
    ALU3dGridLevelIterator(const GridImp & grid, int level, bool);

    //! Constructor for end iterator
    ALU3dGridLevelIterator(const GridImp & grid, int level);

    //! Constructor
    ALU3dGridLevelIterator(const ThisType & org);

    // destructor
    ~ALU3dGridLevelIterator();

    //! prefix increment
    void increment ();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! release entity
    void releaseEntity () {}

    //! assignment of iterators
    ThisType & operator = (const ThisType & org);
  private:
    //! do assignment
    void assign (const ThisType & org);

    // actual level
    int level_;

    // the internal iterator
    IteratorType * iter_ ;

    // deletes iter_
    void removeIter ();

    IteratorType & internalIterator ()
    {
      assert( iter_ );
      return *iter_;
    }
  };

  //********************************************************************
  //
  //  --ALU3dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************
  //! Leaf iterator
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  class ALU3dGridLeafIterator
    : public ALU3dGridEntityPointer<cdim,GridImp>,
      // public LeafIteratorDefaultImplementation<cdim, pitype, GridImp, ALU3dGridLeafIterator>,
      public ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> >
  {
    enum { dim = GridImp :: dimension };

    friend class ALU3dGridEntity<cdim,dim,GridImp>;
    enum { codim = cdim };

    friend class ALU3dGridTreeIterator< ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> > ;

  public:
    typedef typename GridImp::template Codim<cdim>::Entity Entity;

    // the wrapper for the original iterator of the ALU3dGrid
    typedef typename ALU3DSPACE ALU3dGridLeafIteratorWrapper<cdim, pitype> IteratorType;
    typedef IteratorType InternalIteratorType;
    typedef typename ALU3DSPACE IteratorElType<cdim>::val_t val_t;

    typedef ALU3dGridLeafIterator<cdim, pitype, GridImp> ThisType;

    //! Constructor for end iterators
    ALU3dGridLeafIterator(const GridImp & grid, int level);

    //! Constructor for begin Iterators
    ALU3dGridLeafIterator(const GridImp & grid, int level , bool isBegin);

    //! copy Constructor
    ALU3dGridLeafIterator(const ThisType & org);

    //! destructor deleting real iterator
    ~ALU3dGridLeafIterator();

    //! prefix increment
    void increment ();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! release entity
    void releaseEntity () {}

    //! assignment of iterators
    ThisType & operator = (const ThisType & org);

  private:
    // the internal iterator
    IteratorType * iter_;

    // max level for iteration
    int walkLevel_ ;

    //! do assignment
    void assign (const ThisType & org);

    // deletes iter_
    void removeIter () ;

    // return reference to iter_
    InternalIteratorType & internalIterator ()
    {
      assert( iter_ );
      return *iter_;
    }
  };

  // - HierarchicIteraor
  // --HierarchicIterator
  template<class GridImp>
  class ALU3dGridHierarchicIterator
    : public ALU3dGridEntityPointer<0,GridImp>
      // public HierarchicIteratorDefaultImplementation <GridImp,ALU3dGridHierarchicIterator>
  {
    enum { dim = GridImp::dimension };
    typedef ALU3dGridHierarchicIterator<GridImp> ThisType;
    typedef ALU3DSPACE HElementType HElementType ;
    typedef ALU3DSPACE HBndSegType HBndSegType  ;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ctype;

    //! the normal Constructor
    ALU3dGridHierarchicIterator(const GridImp &grid,
                                const HElementType & elem, int maxlevel, bool end );

#ifdef ALU3DGRID_PARALLEL
    //! start constructor for ghosts
    ALU3dGridHierarchicIterator(const GridImp &grid,
                                const HBndSegType& ghost,
                                int maxlevel,
                                bool end);
#endif

    //! the normal Constructor
    ALU3dGridHierarchicIterator(const ALU3dGridHierarchicIterator<GridImp> &org);

    //! increment
    void increment();

    //! dereference Entity, faster then the entity pointersmethod
    Entity & dereference () const;

    //! release entity
    void releaseEntity () {}

    //! the assignment operator
    ThisType & operator = (const ThisType & org);

  private:
    // assign iterator
    void assign(const ThisType & org);

    //! return level of item
    int getLevel(const HElementType* item) const;

    //! return correct level for ghosts
    int getLevel(const HBndSegType* face) const;

    // go to next valid element
    template <class HItemType>
    HItemType* goNextElement (const HItemType* startElem, HItemType * oldEl);

    //! element from where we started
    const HElementType * elem_;

#ifdef ALU3DGRID_PARALLEL
    // pointers to ghost and current ghost
    const HBndSegType * ghost_;
    HBndSegType * nextGhost_;
#endif

    //! maximal level to go down
    int maxlevel_;
  };


} // end namespace Dune

#include "iterator_imp.cc"

#endif
