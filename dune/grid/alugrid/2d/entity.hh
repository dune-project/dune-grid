// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDENTITY_HH
#define DUNE_ALU2DGRIDENTITY_HH

// System includes

// Dune includes
#include <dune/grid/common/entity.hh>
//#include <dune/grid/common/intersectioniteratorwrapper.hh>

// Local includes
#include "iterator.hh"

namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template< int codim, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridGeometry;
  template<class GridImp>
  class ALU2dGridHierarchicIterator;
  template<class GridImp>
  class ALU2dGridLevelIntersectionIterator;
  template<class GridImp>
  class ALU2dGridLeafIntersectionIterator;
  template<class GridImp>
  class ALU2dGridIntersectionIterator;
  template<int codim, PartitionIteratorType, class GridImp>
  class ALU2dGridLeafIterator;
  template<int dim, int dimworld>
  class ALU2dGrid;

  //**********************************************************************
  //
  // --ALU2dGridEntity
  // --Entity
  //
  //**********************************************************************
  /*!
     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

     Here: the general template
   */
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity :
    public EntityDefaultImplementation <cd,dim,GridImp,ALU2dGridEntity>
  {
    static const int dimworld = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    friend class ALU2dGrid < dim , dimworld>;
    friend class ALU2dGridIntersectionIterator < GridImp >;
    friend class ALU2dGridIntersectionIterator < const GridImp >;
    friend class ALU2dGridLevelIntersectionIterator < GridImp >;
    friend class ALU2dGridLevelIntersectionIterator < const GridImp >;
    friend class ALU2dGridLeafIntersectionIterator < GridImp >;
    friend class ALU2dGridLeafIntersectionIterator < const GridImp >;
    friend class ALU2dGridHierarchicIterator   < const GridImp >;
    friend class ALU2dGridHierarchicIterator   < GridImp >;
    friend class ALU2dGridLevelIterator <0,All_Partition,GridImp>;
    friend class ALU2dGridLevelIterator <1,All_Partition,GridImp>;
    friend class ALU2dGridLevelIterator <2,All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <0, All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <1, All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <2, All_Partition,GridImp>;
    friend class ALU2dGridMakeableEntity<0,dim,GridImp>;

    friend class ALU2dGridHierarchicIndexSet<dim,dimworld>;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType ;

  public:
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<cd>::InterfaceType ElementType;
    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<2>::InterfaceType VertexType;

    //! type of our interface entity
    typedef typename GridImp::template Codim<cd>::Entity Entity;
    //! type of corresponding interface geometry
    typedef typename GridImp::template Codim<cd>::Geometry Geometry;
    //! type of geometry implementation
    typedef MakeableInterfaceObject<Geometry> GeometryObj;
    typedef ALU2dGridGeometry<dim-cd,dimworld, const GridImp> GeometryImp;

    //! tpye of EntityPointer
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! level of this element
    int level () const;

    //! Constructor
    ALU2dGridEntity(const GridImp &grid, int level);

    //! Copy Constructor
    ALU2dGridEntity(const ALU2dGridEntity & org);

    //! geometry of this entity
    const Geometry & geometry () const;

    //! return type of geometry
    GeometryType type() const ;

    //! set item pointer to NULL
    void removeElement();

    /*! private methods, but public because of datahandle and template
        arguments of these methods
     */
    //! set element as normal entity
    void setElement(const ElementType &element, int face=-1, int level = -1) const;
    void setElement(const HElementType & el, const VertexType & vx);
    void setElement(const ALU2dGridEntity & org) const {
      setElement(*(org.item_), org.face_);
    }

    //! compare 2 elements
    bool equals ( const ALU2dGridEntity<cd,dim,GridImp> & org ) const;

    //! return partition type of this entity ( see grid.hh )
    //! dummy implementation return InteriorEntity
    PartitionType partitionType() const { return InteriorEntity; }

    /**
       \brief Id of the boundary which is associated with
       the entity, returns 0 for inner entities, arbitrary int otherwise
     */
    int boundaryId () const;

    /*! Location of this vertex within a mesh entity of codimension 0 on the coarse grid.
       This can speed up on-the-fly interpolation for linear conforming elements
       Possibly this is sufficient for all applications we want on-the-fly.
     */
    EntityPointer ownersFather () const;

    //! my position in local coordinates of the owners father
    FieldVector<alu2d_ctype, dim>& positionInOwnersFather () const;

    //! return reference to grid
    const GridImp& grid() const { return grid_; }

    //! return reference to current item
    ElementType& getItem() const
    {
      assert( item_ );
      return *item_;
    }

    // return internal face
    int getFace() const { return face_; }

    //! index is unique within the grid hierachie and per codim
    int getIndex () const;

  private:
    //! returns reference to geometry implementation
    GeometryImp& geoImpl() const {
      return
        grid_.getRealImplementation(geoObj_);
    }

    //! the grid this entity belongs to
    const GridImp &grid_;

    //! corresponding ALU2dGridElement
    mutable ElementType * item_;
    //! the current geometry

    mutable GeometryObj geoObj_;

    mutable int level_;
    mutable int face_;
  };

  /*!
     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

     Entities of codimension 0 ("elements") are defined through template specialization. Note
     that this specialization has an extended interface compared to the general case

     Entities of codimension 0  allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1 with the
     These neighbors are accessed via an iterator. This allows the implementation of
     non-matching meshes. The number of neigbors may be different from the number of faces/edges
     of an element!
   */
  //***********************
  //
  //  --ALU2dGridEntity
  //  --0Entity
  //
  //***********************
  template<int dim, class GridImp>
  class ALU2dGridEntity<0,dim,GridImp>
    : public EntityDefaultImplementation<0,dim,GridImp,ALU2dGridEntity>
  {
    static const int dimworld = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    friend class ALU2dGrid < dim , dimworld>;
    friend class ALU2dGridIntersectionIterator < GridImp >;
    friend class ALU2dGridIntersectionIterator < const GridImp >;
    friend class ALU2dGridLevelIntersectionIterator < GridImp >;
    friend class ALU2dGridLevelIntersectionIterator < const GridImp >;
    friend class ALU2dGridLeafIntersectionIterator < GridImp >;
    friend class ALU2dGridLeafIntersectionIterator < const GridImp >;
    friend class ALU2dGridHierarchicIterator   < const GridImp >;
    friend class ALU2dGridHierarchicIterator   < GridImp >;
    friend class ALU2dGridLevelIterator <0,All_Partition,GridImp>;
    friend class ALU2dGridLevelIterator <1,All_Partition,GridImp>;
    friend class ALU2dGridLevelIterator <2,All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <0, All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <1, All_Partition,GridImp>;
    friend class ALU2dGridLeafIterator <2, All_Partition,GridImp>;
    friend class ALU2dGridMakeableEntity<0,dim,GridImp>;

    friend class ALU2dGridHierarchicIndexSet<dim,dimworld>;

    typedef typename ALU2dImplTraits< dimworld, eltype >::HElementType HElementType ;

  public:
    //! type of our Geometry interface
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    //! type of corresponding interface local geometry
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;
    //! type of our Geometry implementation
    typedef MakeableInterfaceObject<Geometry> GeometryObj;
    typedef ALU2dGridGeometry<dim,dimworld, const GridImp> GeometryImp;

    //! tpye of intersection iterator
    typedef LeafIntersectionIteratorWrapper< GridImp > ALU2dGridLeafIntersectionIteratorType;
    typedef LevelIntersectionIteratorWrapper< GridImp > ALU2dGridLevelIntersectionIteratorType;
    typedef ALU2dGridLeafIntersectionIteratorType ALU2dGridIntersectionIteratorType;

    //! type of entity interface
    typedef typename GridImp::template Codim<0>::Entity Entity;
    //! tpye of entitypointer interface
    //typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef ALU2dGridEntityPointer<0,GridImp> EntityPointer;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    //! Constructor creating empty Entity
    ALU2dGridEntity(const GridImp &grid, int level);

    //! Constructor creating empty Entity
    ALU2dGridEntity(const ALU2dGridEntity & org);

    //! level of this element
    int level () const ;

    //! geometry of this entity
    const Geometry & geometry () const;

    //! return type of geometry
    GeometryType type() const ;

    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
        with codimension cc.
     */
    template<int cc>
    int count () const {
      enum {c = (cc==0) ? 1 : dim+1};
      return c;
    }

    /**
       \brief Id of the boundary which is associated with
       the entity, returns 0 for inner entities, arbitrary int otherwise
     */
    int boundaryId () const {
      // elements are always inside of our Domain
      return 0;
    }

    /*! Intra-level access to intersection with neighboring elements.
       A neighbor is an entity of codimension 0
       which has an entity of codimension 1 in commen with this entity. Access to neighbors
       is provided using iterators. This allows meshes to be nonmatching. Returns iterator
       referencing the first neighbor. */

    // As ibegin() and iend() are deprecated these methods will deliver a LeafIntersectionIterator
    ALU2dGridIntersectionIteratorType ibegin () const
    {
      return ileafbegin();
    }
    //! Reference to one past the last intersection with neighbor
    ALU2dGridIntersectionIteratorType iend () const
    {
      return ileafend();
    }

    ALU2dGridLevelIntersectionIteratorType ilevelbegin () const
    {
      return ALU2dGridLevelIntersectionIteratorType(grid_, *this, this->level(),false);
    }
    ALU2dGridLevelIntersectionIteratorType ilevelend () const
    {
      return ALU2dGridLevelIntersectionIteratorType(grid_, *this, this->level(),true);
    }
    ALU2dGridLeafIntersectionIteratorType ileafbegin () const
    {
      return ALU2dGridLeafIntersectionIteratorType(grid_, *this, this->level(), false);
    }
    ALU2dGridLeafIntersectionIteratorType ileafend () const
    {
      return ALU2dGridLeafIntersectionIteratorType(grid_, *this, this->level(),true);
    }

    //! returns true if Entity is leaf (i.e. has no children)
    bool isLeaf () const;

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    EntityPointer father () const;

    //! returns true if father entity exists
    bool hasFather () const
    {
      return (this->level()>0);
    }

    /*! Inter-level access to son elements on higher levels<=maxlevel.
       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    ALU2dGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
    {
      return ALU2dGridHierarchicIterator<GridImp> (grid_,*item_,maxLevel,false);
    }

    //! Returns iterator to one past the last son
    ALU2dGridHierarchicIterator<GridImp> hend (int maxLevel) const
    {
      return ALU2dGridHierarchicIterator<GridImp> (grid_,*item_,maxLevel,true);
    }

    //! Provide access to mesh entity i of given codimension. Entities
    //!  are numbered 0 ... count<cc>()-1
    template <int cc>
    typename Codim<cc>::EntityPointer entity (int i) const;

    //! Provide access to mesh entity i of given codimension. Entities
    //!  are numbered 0 ... count<cc>()-1
    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( int i ) const
    {
      // apply mapping for codim 1
      return entity< codim >( (codim == 1) ? ( 2 - i ) : i );
    }

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const
    {
#if ALU2DGRID_PARALLEL
      return grid_.rankManager().partitionType( item_->getIndex() );
#else
      return InteriorEntity;
#endif
    }

    /**
       \brief The boundaryId of the i-th subentity of codimension <tt>cc</tt>

       This does the same as <code>entity<cc>(i).boundaryId()</code>, but it is
       usually a lot faster.
     */
    template <int cc>
    int subBoundaryId  ( int i ) const;


    /*! Location of this element relative to the reference element
       of the father. This is sufficient to interpolate all
       dofs in conforming case. Nonconforming may require access to
       neighbors of father and computations with local coordinates.
       On the fly case is somewhat inefficient since dofs  are visited
       several times. If we store interpolation matrices, this is tolerable.
       We assume that on-the-fly implementation of numerical algorithms
       is only done for simple discretizations. Assumes that meshes are nested.
     */

    const LocalGeometry & geometryInFather () const;

    //! The former state() method has been replaced by:
    bool mightVanish () const
    {
      return ((item_->is(ALU2DSPACE Refco::crs))==1);
    }

    bool isNew () const
    {
      return ((item_->wasRefined())==1);
    }

    //***************************************************************
    //  Interface for Adaptation
    //***************************************************************

  public:
    //! marks an element for refCount refines. if refCount is negative the
    //! element is coarsend -refCount times
    //! mark returns true if element was marked, otherwise false
    bool mark(int refCount) const;

    //! return current adaptation mark of element
    int getMark() const;

    /*! private methods, but public because of datahandle and template
        arguments of these methods
     */
    void setElement(const HElementType &element, int face=-1, int level = -1) const;

    void setElement(const ALU2dGridEntity & org) const {
      setElement(*(org.item_));
    }

    //! set actual walk level
    void reset ( int l );

    //! set item pointer to NULL
    void removeElement();

    //! compare 2 entities, which means compare the item pointers
    bool equals ( const ALU2dGridEntity<0,dim,GridImp> & org ) const;

    // return reference to HElement (needed by IntersectionIterator)
    HElementType & getItem() const
    {
      assert( item_ );
      return *item_;
    }

    //! return reference to grid
    const GridImp& grid() const { return grid_; }

    // return internal face
    int getFace() const { return -1; }

    //! index is unique within the grid hierachie and per codim
    int getIndex () const;

  private:
    //! return which number of child we are, i.e. 0 or 1
    int nChild () const;

    //! returns reference to geometry implementation
    GeometryImp& geoImpl() const { return grid_.getRealImplementation(geoObj_); }

    //! return index of sub entity with codim = cc and local number i
    //! i.e. return global number of vertex i
    //! for use in hierarchical index set
    template<int cc> int getSubIndex (int i) const;

    int subIndex (int i, unsigned int codim) const;

    //! corresponding grid
    const GridImp  & grid_;

    //! the current element of grid
    mutable HElementType *item_;

    //! the current geometry
    mutable GeometryObj geoObj_;

    //! is true if entity is leaf entity
    mutable bool isLeaf_;

  }; // end of ALU2dGridEntity codim = 0


  //**********************************************************************
  //
  // --ALU2dGridEntityPointer
  // --EntityPointer
  // --EnPointer
  //**********************************************************************
  /*!
     Enables iteration over all entities of a given codimension and level of a grid.
   */
  template< int codim, class GridImp >
  class ALU2dGridEntityPointer
  {
    // type of this class
    typedef ALU2dGridEntityPointer< codim, GridImp > ThisType;

    static const int dim = GridImp::dimension;
    static const int dimworld = GridImp::dimensionworld;
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    typedef typename Dune::ALU2dImplTraits< dimworld, eltype >::template Codim<codim>::InterfaceType ElementType;

  public:
    enum { codimension = codim };

    //! type of stored entity (interface)
    typedef typename GridImp::template Codim<codimension>::Entity Entity;
    //! tpye of stored entity (implementation)
    typedef ALU2dGridEntity<codimension,dim,GridImp> EntityImp;
    typedef MakeableInterfaceObject<Entity> EntityObj;

    typedef ALU2dGridEntityPointer<codimension,GridImp> EntityPointerImp;

    //! Constructor for EntityPointer that points to an element
    ALU2dGridEntityPointer ( const GridImp &grid,
                             const ElementType &item,
                             int face = -1,
                             int level = -1
                             );

    //! Constructor for EntityPointer init of Level- and LeafIterator
    ALU2dGridEntityPointer(const EntityImp& entity) ;

    //! Constructor for EntityPointer init of Level- and LeafIterator
    ALU2dGridEntityPointer(const GridImp & grid) ;

    //! Copy Constructor
    ALU2dGridEntityPointer(const ThisType & org) ;

    //! Destructor
    ~ALU2dGridEntityPointer();

    //! equality
    // this may have to be changed!
    bool equals (const ThisType & i) const;

    //! dereferencing
    Entity & dereference() const ;

    //! release entity in order to reduce memory allocation
    void compactify();

    //! ask for level of entities
    int level () const;

    ThisType & operator = (const ThisType & org);

  protected:
    EntityImp & entityImp();

    const EntityImp & entityImp() const;

    //! has to be called when iterator is finished
    void done ();

    //! update underlying item pointer and set entity
    void updateEntityPointer( ElementType * item, int face=-1, int level=-1 );

    //! reference to grid
    const GridImp & grid_;

    //! pointer to the real (H)Element
    mutable ElementType * item_;
    mutable int level_;
    int face_;
    //! entity that this EntityPointer points to
    mutable EntityObj * entity_;
  };

} // end namespace Dune

#include "entity_imp.cc"
#endif
