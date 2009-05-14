// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITY_HH
#define DUNE_ALBERTA_ENTITY_HH

#include <dune/grid/common/entity.hh>

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/geometry.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< int codim, class GridImp >
  class AlbertaGridEntityPointer;

  template< int codim, class GridImp, bool leafIterator >
  class AlbertaGridTreeIterator;

  template< class GridImp >
  class AlbertaGridHierarchicIterator;

  template< class GridImp >
  class AlbertaGridIntersection;



  // AlbertaGridEntity
  // -----------------

  /*!
     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

     Here: the general template
   */
  template< int codim, int dim, class GridImp >
  class AlbertaGridEntity
    : public EntityDefaultImplementation< codim, dim, GridImp, AlbertaGridEntity >
  {
    typedef AlbertaGridEntity< codim, dim, GridImp > This;

    enum { dimworld = GridImp::dimensionworld };
    friend class AlbertaGrid< dim , dimworld >;
    friend class AlbertaGridEntity< 0, dim, GridImp>;

    template< int, class, bool > friend class AlbertaGridTreeIterator;
    friend class AlbertaGridEntityPointer< codim, GridImp >;

  public:
    static const int dimension = dim;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    template< int cd >
    struct Codim
    {
      typedef typename GridImp::template Codim< cd >::EntityPointer EntityPointer;
    };

    typedef typename GridImp::template Codim< codim >::Entity Entity;
    typedef typename GridImp::template Codim< codim >::Geometry Geometry;
    typedef typename GridImp::template Codim< codim >::LevelIterator LevelIterator;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  private:
    typedef MakeableInterfaceObject< Geometry > GeometryObject;
    typedef typename GeometryObject::ImplementationType GeometryImp;

  public:
    //! contructor
    AlbertaGridEntity ( const GridImp &grid );

    //! copy constructor
    AlbertaGridEntity ( const This &other );

    //! level of this element
    int level () const;

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const;

    //! geometry of this entity
    const Geometry & geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    //***********************************************
    //  End of Interface methods
    //***********************************************
    //! needed for the LevelIterator and LeafIterator
    ALBERTA EL_INFO *getElInfo () const;
    //! return element for equaltiy in EntityPointer
    ALBERTA EL *getElement () const;

    const ElementInfo &elementInfo () const
    {
      return elementInfo_;
    }

    //! equality of entities
    bool equals ( const This &other ) const;

    void clearElement ();
    void setElement ( const ElementInfo &elementInfo, int subEntity );

    // same as setElInfo just with a entity given
    void setEntity ( const This &other );

    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return grid_;
    }

    //! obtain number of the subentity within the element (in ALBERTA numbering)
    int subEntity () const
    {
      return subEntity_;
    }

    //! obtain twist
    int twist () const
    {
      return elementInfo().template twist< codimension >( subEntity() );
    }

  private:
    const GeometryImp &geoImp () const
    {
      return GridImp :: getRealImplementation( geo_ );
    }

    GeometryImp &geoImp ()
    {
      return GridImp :: getRealImplementation( geo_ );
    }

  private:
    // grid this entity belong to
    const GridImp &grid_;

    // ALBERTA element info
    ElementInfo elementInfo_;

    // number of the subentity within the element (in ALBERTA numbering)
    int subEntity_;

    // current geometry
    GeometryObject geo_;
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
  //  --AlbertaGridEntity
  //  --0Entity
  //
  //***********************
  template< int dim, class GridImp >
  class AlbertaGridEntity< 0, dim, GridImp >
    : public EntityDefaultImplementation< 0, dim, GridImp, AlbertaGridEntity >
  {
    typedef AlbertaGridEntity< 0, dim, GridImp > This;

    static const int dimworld = GridImp::dimensionworld;

    friend class AlbertaGrid< dim, dimworld >;
    friend class AlbertaGridIntersection< GridImp >;
    friend class AlbertaGridHierarchicIterator< GridImp >;
    template< int, class, bool > friend class AlbertaGridTreeIterator;
    friend class AlbertaGridEntityPointer<0,GridImp>;

  public:
    static const int dimension = dim;
    static const int codimension = 0;
    static const int mydimension = dimension - codimension;

    template< int codim >
    struct Codim
    {
      typedef typename GridImp::template Codim< codim >::EntityPointer
      EntityPointer;
    };

    typedef typename GridImp::template Codim< 0 >::Entity Entity;
    typedef typename GridImp::template Codim< 0 >::Geometry Geometry;
    typedef typename GridImp::template Codim< 0 >::LocalGeometry LocalGeometry;
    typedef AlbertaGridGeometry< dimension, dimworld, GridImp > GeometryImp;

    typedef typename GridImp::template Codim<0>::LevelIterator LevelIterator;
    typedef typename GridImp::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    typedef AlbertaGridLeafIntersectionIterator< GridImp > AlbertaGridLeafIntersectionIteratorType;
    typedef AlbertaGridLeafIntersectionIteratorType AlbertaGridIntersectionIteratorType;
    typedef AlbertaGridLeafIntersectionIteratorType AlbertaGridLevelIntersectionIteratorType;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

    //! Destructor, needed perhaps needed for deleteing faceEntity_ and
    //! edgeEntity_ , see below
    //! there are only implementations for dim==dimworld 2,3
    ~AlbertaGridEntity() {};

    //! Constructor, real information is set via setElInfo method
    explicit AlbertaGridEntity ( const GridImp &grid );

    AlbertaGridEntity ( const This &other );

    //! level of this element
    int level () const;

    //! index of the boundary which is associated with the entity, 0 for inner entities
    int boundaryId () const;

    //! geometry of this entity
    const Geometry & geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    /** obtain the number of subentities of a codimension
     *
     *  \tparam  codim  codimension
     *
     *  \returns the number of subentities of the given codimension
     */
    template< int codim >
    int count () const
    {
      return Alberta::NumSubEntities< dimension, codim >::value;
    }

    /** obtain a subentity
     *
     *  \tparam  codim  codimension of the desired subentity
     *
     *  \param[in]  i  number of the subentity (in generic numbering)
     *
     *  \returns an EntityPointer to the subentity
     *
     *  \note: The subentities are numbered 0, ..., count< codim >-1
     */
    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( int i ) const;

    /*! Intra-level access to intersection with neighboring elements.
       A neighbor is an entity of codimension 0
       which has an entity of codimension 1 in commen with this entity. Access to neighbors
       is provided using iterators. This allows meshes to be nonmatching. Returns iterator
       referencing the first neighbor. */
    AlbertaGridLeafIntersectionIteratorType ileafbegin () const;

    //! Reference to one past the last intersection with neighbor
    AlbertaGridIntersectionIteratorType ileafend () const;

    AlbertaGridLevelIntersectionIteratorType ilevelbegin () const
    {
      DUNE_THROW( NotImplemented, "method ilevelbegin not implemented for AlbertaGrid." );
      return ileafend();
    }

    AlbertaGridLeafIntersectionIteratorType ilevelend () const
    {
      DUNE_THROW( NotImplemented, "method ilevelend not implemented for AlbertaGrid." );
      return ileafend();
    }

    //! returns true if entity is leaf entity, i.e. has no children
    bool isLeaf () const;

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    EntityPointer father () const;

    /** \brief Location of this element relative to the father's reference element
     *
     *  This information is sufficient to interpolate all dofs in conforming case.
     *  Nonconforming may require access to neighbors of father and computations
     *  with local coordinates.
     *  On the fly case is somewhat inefficient since dofs  are visited several
     *  times. If we store interpolation matrices, this is tolerable.
     */
    const LocalGeometry &geometryInFather () const;

    /*! Inter-level access to son elements on higher levels<=maxlevel.
       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    HierarchicIterator hbegin (int maxlevel) const;

    //! Returns iterator to one past the last son
    HierarchicIterator hend (int maxlevel) const;

    /** \brief Was the entity created during the last adaptation cycle? */
    bool isNew () const;

    /**\brief Might the entity vanish during the next adaptation cycle? */
    bool mightVanish () const;

    /**\brief Returns true, if entity has intersections with boundary
     */
    bool hasBoundaryIntersections () const ;

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const;

    //! equality of entities
    bool equals ( const AlbertaGridEntity<0,dim,GridImp> & i) const;

    //***************************************************************
    //  Interface for parallelisation
    //***************************************************************
    // set leaf data with processor number
    void setLeafData( int proc );

    // return true if this entity belong to master set of this grid
    bool master() const;

    // needed for LevelIterator to compare
    ALBERTA EL_INFO *getElInfo () const;

    const ElementInfo &elementInfo () const
    {
      return elementInfo_;
    }

    // return element for equaltiy in EntityPointer
    ALBERTA EL *getElement () const;

    void clearElement ();
    void setElement ( const ElementInfo &elementInfo, int subEntity );

    // same as setElInfo just with a entity given
    void setEntity ( const This &other);

    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return grid_;
    }

    //! obtain number of the subentity within the element (in ALBERTA numbering)
    int subEntity () const
    {
      return 0;
    }

    //! obtain twist
    int twist () const
    {
      return elementInfo().template twist< codimension >( subEntity() );
    }

    //! obtain twist of a subentity
    template< int codim >
    int twist ( int i )
    {
      return elementInfo().template twist< codim >( grid_.dune2alberta( codim, i ) );
    }

  private:
    //! return which number of child we are, i.e. 0 or 1
    int nChild () const;

    GeometryImp &geoImp () const
    {
      return GridImp::getRealImplementation( geo_ );
    }

    //! the corresponding grid
    const GridImp & grid_;

    // Alberta element info
    ElementInfo elementInfo_;

    // local coordinates within father
    typedef MakeableInterfaceObject< Geometry > GeometryObject;

    //! the cuurent geometry
    mutable GeometryObject geo_;
    mutable bool builtgeometry_;  //!< true if geometry has been constructed
  }; // end of AlbertaGridEntity codim = 0

}

#endif
