// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDENTITY_HH
#define DUNE_ALU3DGRIDENTITY_HH

// System includes

// Dune includes
#include <dune/grid/common/entity.hh>
#include <dune/grid/alugrid/common/intersectioniteratorwrapper.hh>

// Local includes
#include "alu3dinclude.hh"
#include "iterator.hh"
#include "entityseed.hh"

namespace Dune
{

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
  template<int codim, PartitionIteratorType, class GridImp>
  class ALU3dGridLeafIterator;
  template< ALU3dGridElementType, class >
  class ALU3dGrid;

  /*!
     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

     Here: the general template
   */
  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity :
    public EntityDefaultImplementation <cd,dim,GridImp,ALU3dGridEntity>
  {
    // default just returns level
    template <class GridType, int cdim>
    struct GetLevel
    {
      template <class ItemType>
      static int getLevel(const GridType & grid, const ItemType & item )
      {
        return item.level();
      }
    };

    // for leaf vertices the level is somewhat difficult to obtain, because
    // this the maximum of levels of elements that have this vertex as sub
    // entity
    template <class GridType>
    struct GetLevel<GridType,3>
    {
      template <class ItemType>
      static int getLevel(const GridType & grid, const ItemType & item)
      {
        return (item.isLeafEntity()) ? grid.getLevelOfLeafVertex(item) : item.level();
      }
    };

    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGrid< GridImp::elementType, Comm >;
    friend class ALU3dGridEntity < 0, dim, GridImp >;
    friend class ALU3dGridLevelIterator < cd, All_Partition, GridImp >;

    friend class ALU3dGridHierarchicIndexSet< GridImp::elementType, Comm >;

    template< class > friend class ALU3dGridFactory;

  public:
    typedef typename GridImp::GridObjectFactoryType FactoryType;

    typedef ALU3dImplTraits< GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<cd>::InterfaceType HItemType;
    typedef typename ImplTraits::template Codim<cd>::ImplementationType ItemType;
    typedef typename ImplTraits::VertexType VertexType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

    typedef typename GridImp::template Codim<cd>::Entity Entity;
    typedef typename GridImp::template Codim<cd>::Geometry Geometry;

    typedef ALU3dGridGeometry<dim-cd,GridImp::dimensionworld,GridImp> GeometryImp;
    typedef MakeableInterfaceObject<Geometry> GeometryObject;

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! typedef of my type
    typedef typename GridImp::template Codim<cd>::EntitySeed EntitySeed;

    //! level of this element
    int level () const;

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const;

    //! Constructor
    ALU3dGridEntity(const FactoryType &factory, int level);

    //! copy Constructor
    ALU3dGridEntity(const ALU3dGridEntity & org);

    //! geometry of this entity
    const Geometry & geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    // set element as normal entity
    void setElement(const HItemType & item);
    void setElement(const HItemType & item, const int level, int twist=0, int face = -1);

    /* set entity from seed */
    void setElement(const EntitySeed& seed);

    //! setGhost is not valid for this codim
    void setGhost(const HBndSegType  &ghost);

    //! reset item pointer to NULL
    void removeElement ();

    //! reset item pointer to NULL
    void reset ( int l );

    //! compare 2 elements by comparing the item pointers
    bool equals ( const ALU3dGridEntity<cd,dim,GridImp> & org ) const;

    //! set item from other entity, mainly for copy constructor of entity pointer
    void setEntity ( const ALU3dGridEntity<cd,dim,GridImp> & org );

    // return reference to internal item
    const ItemType & getItem () const { return *item_; }

    //! return reference to grid
    const GridImp& grid() const { return factory_.grid(); }

    //! return reference to factory
    const FactoryType& factory() const { return factory_; }

    //! return seed of entity
    EntitySeed seed() const
    {
      return EntitySeed( getItem(), level(), twist_, face_ );
    }

  private:
    //! return reference to geometry implementation
    GeometryImp& geoImp() const { return GridImp :: getRealImplementation(geo_); }

    //! index is unique within the grid hierarchy and per codim
    int getIndex () const;

    //! convert ALUGrid partition type to dune partition type
    PartitionType convertBndId(const HItemType & item) const ;

    //! the cuurent geometry
    mutable GeometryObject geo_;

    // the factory that created this entity
    const FactoryType& factory_;

    // corresponding ALU3dGrid item, here face, edge, or vertex
    const ItemType * item_;

    int level_; //! level of entity
    int gIndex_; //! hierarchic index
    int twist_; //! twist of the underlying ALU element (with regard to the element that asked for it)
    int face_; //! for face, know on which face we are

    mutable bool builtgeometry_;        //!< true if geometry has been constructed
    mutable PartitionType partitionType_; //!< partition type of this entity
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
  //  --ALU3dGridEntity
  //  --0Entity
  //
  //***********************
  template<int dim, class GridImp>
  class ALU3dGridEntity<0,dim,GridImp>
    : public EntityDefaultImplementation<0,dim,GridImp,ALU3dGridEntity>
  {
    static const int dimworld = remove_const< GridImp >::type::dimensionworld;
    static const ALU3dGridElementType elementType = remove_const< GridImp >::type::elementType;

    typedef typename GridImp::MPICommunicatorType Comm;

    typedef ALU3dImplTraits< elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<0>::InterfaceType HElementType;

    typedef typename ImplTraits::GEOElementType GEOElementType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

    enum { refine_element_t = ImplTraits::refine_element_t };
    enum { coarse_element_t = ImplTraits::coarse_element_t };
    enum { nosplit_element_t = ImplTraits::nosplit_element_t };

    typedef typename ImplTraits::MarkRuleType MarkRuleType;

    friend class ALU3dGrid< elementType, Comm >;
    friend class ALU3dGridIntersectionIterator < GridImp >;
    friend class ALU3dGridIntersectionIterator < const GridImp >;
    friend class ALU3dGridHierarchicIterator   < const GridImp >;
    friend class ALU3dGridHierarchicIterator   < GridImp >;
    friend class ALU3dGridLevelIterator <0,All_Partition,GridImp>;
    friend class ALU3dGridLevelIterator <1,All_Partition,GridImp>;
    friend class ALU3dGridLevelIterator <2,All_Partition,GridImp>;
    friend class ALU3dGridLevelIterator <3,All_Partition,GridImp>;
    friend class ALU3dGridLeafIterator <0, All_Partition,GridImp>;
    friend class ALU3dGridLeafIterator <1, All_Partition,GridImp>;
    friend class ALU3dGridLeafIterator <2, All_Partition,GridImp>;
    friend class ALU3dGridLeafIterator <3, All_Partition,GridImp>;

    friend class ALU3dGridHierarchicIndexSet< elementType, Comm >;

    template< class > friend class ALU3dGridFactory;

    // type of reference element
    typedef typename GridImp :: ReferenceElementType ReferenceElementType;

  public:
    typedef typename GridImp::GridObjectFactoryType FactoryType;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    typedef ALU3dGridGeometry<dim,dimworld,GridImp> GeometryImp;
    typedef MakeableInterfaceObject<Geometry> GeometryObject;
    typedef ALU3dGridIntersectionIterator<GridImp> IntersectionIteratorImp;

    typedef LeafIntersectionIteratorWrapper <GridImp>  ALU3dGridIntersectionIteratorType;
    typedef LeafIntersectionIteratorWrapper <GridImp>  ALU3dGridLeafIntersectionIteratorType;
    typedef LevelIntersectionIteratorWrapper<GridImp>  ALU3dGridLevelIntersectionIteratorType;

    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    //! typedef of my type
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    //! Constructor creating empty Entity
    ALU3dGridEntity(const FactoryType& factory, int level);

    //! copy Constructor
    ALU3dGridEntity(const ALU3dGridEntity & org);

    //! level of this element
    int level () const ;

    //! geometry of this entity
    const Geometry & geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const;

    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
        with codimension cc.
     */
    template<int cc> int count () const ;

    //! Provide access to mesh entity i of given codimension. Entities
    //!  are numbered 0 ... count<cc>()-1
    template <int codim>
    typename Codim< codim >::EntityPointer entity (int i) const
    {
      typedef GenericGeometry::MapNumberingProvider< GridImp::dimension > Numbering;
      const unsigned int tid = type().id();
      const int j = Numbering::template dune2generic< codim >( tid, i );
      return subEntity< codim >( j );
    }

    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( int i ) const;

    /*! Access to intersection with neighboring elements that are leaf
     * elements. A neighbor is an entity of codimension 0
       which has an entity of codimension 1 in commen with this entity. Access to neighbors
       is provided using iterators. This allows meshes to be nonmatching. Returns iterator
       referencing the first neighbor. */
    ALU3dGridLeafIntersectionIteratorType ileafbegin () const;

    //! Reference to one past the last intersection with neighbor
    ALU3dGridLeafIntersectionIteratorType ileafend () const;

    /*! Intra-level access to intersection with neighboring elements.
       A neighbor is an entity of codimension 0
       which has an entity of codimension 1 in commen with this entity. Access to neighbors
       is provided using iterators. This allows meshes to be nonmatching. Returns iterator
       referencing the first neighbor. */
    ALU3dGridLevelIntersectionIteratorType ilevelbegin () const;

    //! Reference to one past the last intersection with neighbor
    ALU3dGridLevelIntersectionIteratorType ilevelend () const;

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

    /*! Location of this element relative to the reference element
       of the father. This is sufficient to interpolate all
       dofs in conforming case. Nonconforming may require access to
       neighbors of father and computations with local coordinates.
       On the fly case is somewhat inefficient since dofs  are visited
       several times. If we store interpolation matrices, this is tolerable.
       We assume that on-the-fly implementation of numerical algorithms
       is only done for simple discretizations. Assumes that meshes are nested.
     */
    const Geometry & geometryInFather () const;

    /*! Inter-level access to son elements on higher levels<=maxlevel.
       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    ALU3dGridHierarchicIterator<GridImp> hbegin (int maxlevel) const;

    //! Returns iterator to one past the last son
    ALU3dGridHierarchicIterator<GridImp> hend (int maxlevel) const;

    //***************************************************************
    //  Interface for Adaptation
    //***************************************************************

    //! returns true, if entity was created during last adaptation cycle
    bool isNew () const;

    //! returns true, if entity might be coarsened during next adaptation cycle
    bool mightVanish () const;

    //! returns true, if entity has intersections with boundary
    bool hasBoundaryIntersections () const;

    // private method
    //! marks an element for refCount refines. if refCount is negative the
    //! element is coarsend -refCount times
    //! mark returns true if element was marked, otherwise false
    bool mark( int refCount ) const;

    //! \brief return current adaptation mark for this entity
    int getMark() const;

    /*! private methods, but public because of datahandle and template
        arguments of these methods
     */
    void setElement(HElementType &element);

    /* set entity from seed */
    void setElement(const EntitySeed& seed);

    //! set original element pointer to fake entity
    void setGhost(HBndSegType & ghost);

    //! set actual walk level
    void reset ( int l );

    //! set item pointer to NULL
    void removeElement();

    //! compare 2 entities, which means compare the item pointers
    bool equals ( const ALU3dGridEntity<0,dim,GridImp> & org ) const;

    void setEntity ( const ALU3dGridEntity<0,dim,GridImp> & org );

    //! return index of sub entity with codim = cc and local number i
    //! i.e. return global number of vertex i
    //! for use in hierarchical index set
    template<int cc> int getSubIndex (int i) const;

    //! return index of sub entity with codim = cc and local number i
    //! i.e. return global number of vertex i
    //! for use in hierarchical index set
    int subIndex(int i, unsigned int codim) const;

    // return reference to internal item
    const IMPLElementType& getItem () const { return *item_; }

    // return reference to internal item
    const BNDFaceType& getGhost () const
    {
      assert( isGhost() );
      return *ghost_;
    }

    //! return reference to grid
    const GridImp& grid() const { return factory_.grid(); }

    //! return reference to factory
    const FactoryType& factory() const { return factory_; }

    //! returns true if entity is ghost
    bool isGhost () const { return ImplTraits::isGhost( ghost_ ); }

    //! return key for this entity
    EntitySeed seed() const
    {
      if( isGhost() )
        return EntitySeed( getGhost () );
      else
        return EntitySeed( getItem() );
    }

  private:
    //! return reference to geometry implementation
    GeometryImp& geoImp() const { return GridImp :: getRealImplementation(geo_); }

    //! index is unique within the grid hierachy and per codim
    int getIndex () const;

    //! the entity's geometry
    mutable GeometryObject geo_;

    // corresponding factory
    const FactoryType& factory_;

    // the current element of grid
    mutable IMPLElementType* item_;

    //! not zero if entity is ghost entity
    mutable BNDFaceType*  ghost_;

    int level_;  //!< level of element
    bool isLeaf_; //!< is true if entity is leaf entity
    mutable bool builtgeometry_; //!< true if geometry has been constructed

  }; // end of ALU3dGridEntity codim = 0
     //**********************************************************************
     //
     // --ALU3dGridEntityPointer
     // --EntityPointer
     // --EnPointer
     /*!
        Enables iteration over all entities of a given codimension and level of a grid.
      */
  template< int codim, class GridImp >
  class ALU3dGridEntityPointerBase
  //: public EntityPointerDefaultImplementation <codim, GridImp, ALU3dGridEntityPointer<cd,GridImp> >
  {
    typedef ALU3dGridEntityPointerBase< codim, GridImp > ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity<codim,dim,GridImp>;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits<GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<codim>::InterfaceType HElementType;

    typedef typename ImplTraits::HBndSegType HBndSegType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;
  public:
    typedef typename GridImp::GridObjectFactoryType FactoryType;

    enum { codimension = codim };

    //! type of Entity
    typedef typename GridImp::template Codim<codimension>::Entity Entity;
    //! underlying EntityImplementation
    typedef MakeableInterfaceObject<Entity> EntityObject;
    typedef typename EntityObject :: ImplementationType EntityImp;

    //! typedef of my type
    typedef ThisType ALU3dGridEntityPointerType;

    //! make type of entity pointer implementation available in derived classes
    typedef ALU3dGridEntityPointer<codimension,GridImp> EntityPointerImp;

    //! type of entity seed
    typedef ALU3dGridEntitySeed<codimension, GridImp> ALU3dGridEntitySeedType;

    //! Constructor for EntityPointer that points to an element
    ALU3dGridEntityPointerBase(const FactoryType& factory,
                               const HElementType & item);

    //! Constructor for EntityPointer that points to an element
    ALU3dGridEntityPointerBase(const FactoryType& factory,
                               const HElementType & item,
                               const int level,
                               const int twist,
                               const int duneFace
                               );

    //! Constructor for EntityPointer that points to an ghost
    ALU3dGridEntityPointerBase(const FactoryType& factory,
                               const HBndSegType & ghostFace );

    //! Constructor for EntityPointer that points to an ghost
    ALU3dGridEntityPointerBase(const FactoryType& factory,
                               const ALU3dGridEntitySeedType& seed );

    //! copy constructor
    ALU3dGridEntityPointerBase(const ALU3dGridEntityPointerType & org);

    //! Destructor
    ~ALU3dGridEntityPointerBase();

    //! equality
    bool equals (const ALU3dGridEntityPointerType& i) const;

    //! assignment operator
    ThisType & operator = (const ThisType & org);

    //! dereferencing
    Entity & dereference () const ;

    //! ask for level of entities
    int level () const ;

  protected:
    // clones object
    void clone (const ALU3dGridEntityPointerType & org);

    // get entity and assign from org.entity
    void getEntity (const ALU3dGridEntityPointerType & org);

    //! has to be called when iterator is finished
    void done ();

    //! put entity to entity stack
    void freeEntity ();

    //! return reference to grid
    const GridImp& grid () const { return factory_.grid(); }

    //! Constructor for EntityPointer init of Level-, and Leaf-, and
    //! HierarchicIterator
    ALU3dGridEntityPointerBase(const FactoryType& factory, int level );

    // update underlying item pointer and set ghost entity
    void updateGhostPointer( HBndSegType & ghostFace );
    // update underlying item pointer and set entity
    void updateEntityPointer( HElementType * item , int level = -1 );

    // reference to factory
    const FactoryType& factory_;

    // key to gererate entity
    ALU3dGridEntitySeedType seed_;

    // entity that this EntityPointer points to
    mutable EntityObject * entity_;

    // is true if entity must not be released
    bool locked_;

    // return reference to internal entity implementation
    EntityImp & entityImp () const {
      assert( entity_ );
      return GridImp :: getRealImplementation(*entity_);
    }
  };

  //! ALUGridEntityPointer points to an entity
  //! this class is the specialisation for codim 0,
  //! it has exactly the same functionality as the ALU3dGridEntityPointerBase
  template<class GridImp>
  class ALU3dGridEntityPointer<0,GridImp> :
    public ALU3dGridEntityPointerBase<0,GridImp>
  {
  protected:
    typedef ALU3dGridEntityPointerBase<0,GridImp> BaseType;

    enum { cd = 0 };
    typedef ALU3dGridEntityPointer <cd,GridImp> ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity<cd,dim,GridImp>;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits<GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<cd>::InterfaceType HElementType;

    typedef typename ImplTraits::HBndSegType HBndSegType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

    typedef ALU3dGridEntity< 0,dim,GridImp> ALU3dGridEntityType ;

    using BaseType :: seed_;
    using BaseType :: entity_;
    using BaseType :: entityImp;
    using BaseType :: factory_;
  public:
    typedef typename GridImp::GridObjectFactoryType FactoryType;

    //! type of entity seed
    typedef ALU3dGridEntitySeed<cd, GridImp> ALU3dGridEntitySeedType;

    //! type of Entity
    typedef typename GridImp::template Codim<cd>::Entity Entity;

    //! typedef of my type
    typedef ThisType ALU3dGridEntityPointerType;

    //! Constructor for EntityPointer that points to an interior element
    ALU3dGridEntityPointer(const FactoryType& factory,
                           const HElementType & item)
      : ALU3dGridEntityPointerBase<cd,GridImp> (factory,item) {}

    //! Constructor for EntityPointer that points to an ghost
    ALU3dGridEntityPointer(const FactoryType& factory,
                           const HBndSegType & ghostFace )
      : ALU3dGridEntityPointerBase<cd,GridImp> (factory,ghostFace) {}

    //! Constructor for EntityPointer that points to given entity
    ALU3dGridEntityPointer(const FactoryType& factory, const ALU3dGridEntitySeedType& seed)
      : ALU3dGridEntityPointerBase<cd,GridImp> (factory, seed)
    {
      // for ghost entities we have to copy right away
      if( seed.isGhost() )
      {
        assert( entity_ == 0 );
        entity_ = factory_.template getNewEntity<0> ();
        assert( entity_ );
        entityImp().setGhost( *seed.ghost() );

        // don't free on compactify, otherwise ghost info is lost
        this->locked_ = true ;
      }
    }

    //! Constructor for EntityPointer that points to an entity (interior or ghost)
    ALU3dGridEntityPointer(const ALU3dGridEntityType& entity)
      : ALU3dGridEntityPointerBase<cd,GridImp> (entity.factory(),
                                                entity.seed() )
    {
      // for ghost entities we have to copy right away
      if( entity.isGhost() )
      {
        assert( entity_ == 0 );
        entity_ = factory_.template getNewEntity<0> ();
        assert( entity_ );
        entityImp().setEntity( entity );

        // don't free on compactify, otherwise ghost info is lost
        this->locked_ = true ;
      }
    }

    //! copy constructor
    ALU3dGridEntityPointer(const ALU3dGridEntityPointerType & org)
      : ALU3dGridEntityPointerBase<cd,GridImp> (org)
    {}

  protected:
    //! Constructor for EntityPointer init of Level-, and Leaf-, and
    //! HierarchicIterator
    ALU3dGridEntityPointer(const FactoryType& factory, int level )
      : ALU3dGridEntityPointerBase<cd,GridImp> (factory,level) {}
  };


  template<int cd, class GridImp>
  class ALU3dGridEntityPointer :
    public ALU3dGridEntityPointerBase<cd,GridImp>
  {
  protected:
    typedef ALU3dGridEntityPointerBase<cd,GridImp> BaseType ;
    typedef ALU3dGridEntityPointer <cd,GridImp> ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity<cd,dim,GridImp>;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits< GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<cd>::InterfaceType HElementType;

    typedef typename ImplTraits::HBndSegType HBndSegType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;
    typedef ALU3dGridEntity<cd,dim,GridImp> ALU3dGridEntityType;

    using BaseType :: seed_;
    using BaseType :: entity_;
    using BaseType :: entityImp;
    using BaseType :: factory_;
    using BaseType :: getEntity;

  public:
    typedef typename GridImp::GridObjectFactoryType FactoryType;

    //! type of entity seed
    typedef ALU3dGridEntitySeed<cd, GridImp> ALU3dGridEntitySeedType;

    //! type of Entity
    typedef typename GridImp::template Codim<cd>::Entity Entity;

    //! typedef of my type
    typedef ALU3dGridEntityPointer<cd,GridImp> ALU3dGridEntityPointerType;

  protected:
    static const int defaultValue = -665; //ALU3dGridEntityPointerType :: defaultValue;

  public:
    //! Constructor for EntityPointer that points to an element
    ALU3dGridEntityPointer(const FactoryType& factory,
                           const int level,
                           const HElementType & item,
                           const int twist = defaultValue,
                           const int duneFace = defaultValue
                           );

    //! Constructor for EntityPointer that points to given entity
    ALU3dGridEntityPointer(const ALU3dGridEntityType& entity)
      : ALU3dGridEntityPointerBase<cd,GridImp> (entity.factory(),
                                                entity.seed())
    {}

    //! Constructor for EntityPointer that points to given entity
    ALU3dGridEntityPointer(const FactoryType& factory, const ALU3dGridEntitySeedType& seed)
      : ALU3dGridEntityPointerBase<cd,GridImp> ( factory, seed )
    {}

    //! copy constructor
    ALU3dGridEntityPointer(const ALU3dGridEntityPointerType & org);

    //! dereferencing
    Entity & dereference () const ;

    //! assignment operator
    ThisType & operator = (const ThisType & org);

    //! ask for level of entities
    int level () const ;

  protected:
    // clones object
    void clone (const ALU3dGridEntityPointerType & org);

    void updateEntityPointer( HElementType * item , int level );

    //! Constructor for EntityPointer init of Level-, and Leaf-, and
    //! HierarchicIterator
    ALU3dGridEntityPointer(const FactoryType& factory, int level )
      : ALU3dGridEntityPointerBase<cd,GridImp> (factory,level)
    {}
  };

} // end namespace Dune

#include "entity_inline.hh"

#if COMPILE_ALUGRID_INLINE
  #include "entity_imp.cc"
#endif

#endif
