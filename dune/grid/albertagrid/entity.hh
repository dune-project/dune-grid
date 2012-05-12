// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITY_HH
#define DUNE_ALBERTA_ENTITY_HH

#include <dune/grid/common/entity.hh>

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/entityseed.hh>
#include <dune/grid/albertagrid/geometry.hh>

#if HAVE_ALBERTA

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
  class AlbertaGridLeafIntersection;

  template< class GridImp >
  class AlbertaGridLeafIntersectionIterator;



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
    typedef typename GridImp::template Codim< codim >::EntitySeed EntitySeed;
    typedef typename GridImp::template Codim< codim >::Geometry Geometry;
    typedef typename GridImp::template Codim< codim >::LevelIterator LevelIterator;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  private:
    //typedef MakeableInterfaceObject< Geometry > GeometryObject;
    //typedef typename GeometryObject::ImplementationType GeometryImp;
    typedef typename GridImp::Traits::template Codim< codim >::GeometryImpl GeometryImpl;

  public:
    //! constructor
    explicit AlbertaGridEntity ( const GridImp &grid );

    //! contructor
    AlbertaGridEntity ( const GridImp &grid, const ElementInfo &elementInfo, int subEntity );

    //! level of this element
    int level () const;

    //! return partition type of this entity ( see grid.hh )
    PartitionType partitionType() const;

    //! geometry of this entity
    Geometry geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    //! obtain entity seed
    EntitySeed seed () const { return EntitySeed( elementInfo(), subEntity() ); }

    //***********************************************
    // end of interface methods
    //***********************************************

    //! needed for the LevelIterator and LeafIterator
    ALBERTA EL_INFO *getElInfo () const;

    const ElementInfo &elementInfo () const { return elementInfo_; }

    //! equality of entities
    bool equals ( const This &other ) const;

    void clearElement ();
    void setElement ( const ElementInfo &elementInfo, int subEntity );

    // same as setElInfo just with a entity given
    void setEntity ( const This &other );

    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return *grid_;
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
    // grid this entity belong to
    const GridImp *grid_;

    // ALBERTA element info
    ElementInfo elementInfo_;

    // number of the subentity within the element (in ALBERTA numbering)
    int subEntity_;

    // current geometry
    GeometryImpl geo_;
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
  template< int dim, class GridImp >
  class AlbertaGridEntity< 0, dim, GridImp >
    : public EntityDefaultImplementation< 0, dim, GridImp, AlbertaGridEntity >
  {
    typedef AlbertaGridEntity< 0, dim, GridImp > This;

    static const int dimworld = GridImp::dimensionworld;

    friend class AlbertaGrid< dim, dimworld >;
    friend class AlbertaGridLeafIntersection< GridImp >;
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
    typedef typename GridImp::template Codim< 0 >::EntitySeed EntitySeed;
    typedef typename GridImp::template Codim< 0 >::Geometry Geometry;
    typedef typename GridImp::template Codim< 0 >::LocalGeometry LocalGeometry;
    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;

    typedef typename GridImp::template Codim<0>::LevelIterator LevelIterator;
    typedef typename GridImp::HierarchicIterator HierarchicIterator;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    typedef Dune::AlbertaGridLeafIntersectionIterator< GridImp > AlbertaGridLeafIntersectionIterator;
    typedef AlbertaGridLeafIntersectionIterator AlbertaGridLevelIntersectionIterator;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

    //! constructor
    explicit AlbertaGridEntity ( const GridImp &grid );

    //! constructor
    AlbertaGridEntity ( const GridImp &grid, const ElementInfo &elementInfo, int subEntity );

    //! level of this element
    int level () const;

    //! index of the boundary which is associated with the entity, 0 for inner entities
    int boundaryId () const;

    //! geometry of this entity
    Geometry geometry () const;

    //! type of geometry of this entity
    GeometryType type () const;

    //! obtain entity seed
    EntitySeed seed () const { return EntitySeed( elementInfo() ); }

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
    AlbertaGridLeafIntersectionIterator ileafbegin () const;

    //! Reference to one past the last intersection with neighbor
    AlbertaGridLeafIntersectionIterator ileafend () const;

    AlbertaGridLevelIntersectionIterator ilevelbegin () const
    {
      if( grid().maxLevel() == 0 )
        return ileafbegin();
      else
      {
        DUNE_THROW( NotImplemented, "method ilevelbegin not implemented for AlbertaGrid." );
        return ileafend();
      }
    }

    AlbertaGridLevelIntersectionIterator ilevelend () const
    {
      return ileafend();
    }

    //! returns true if entity is leaf entity, i.e. has no children
    bool isLeaf () const;

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    EntityPointer father () const;
    //! returns true if father entity exists
    bool hasFather () const
    {
      return (this->level()>0);
    }

    /** \brief Location of this element relative to the father's reference element
     *
     *  This information is sufficient to interpolate all dofs in conforming case.
     *  Nonconforming may require access to neighbors of father and computations
     *  with local coordinates.
     *  On the fly case is somewhat inefficient since dofs  are visited several
     *  times. If we store interpolation matrices, this is tolerable.
     */
    LocalGeometry geometryInFather () const;

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

    // needed for LevelIterator to compare
    ALBERTA EL_INFO *getElInfo () const;

    const ElementInfo &elementInfo () const
    {
      return elementInfo_;
    }

    void clearElement ();
    void setElement ( const ElementInfo &elementInfo, int subEntity );

    // same as setElInfo just with a entity given
    void setEntity ( const This &other );

    //! obtain a reference to the grid
    const GridImp &grid () const
    {
      return *grid_;
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
    int twist ( int i ) const
    {
      return elementInfo().template twist< codim >( grid().generic2alberta( codim, i ) );
    }

  private:
    //! return which number of child we are, i.e. 0 or 1
    int nChild () const;

    //! the corresponding grid
    const GridImp *grid_;

    // Alberta element info
    ElementInfo elementInfo_;

    // local coordinates within father
    typedef MakeableInterfaceObject< Geometry > GeometryObject;

    //! the cuurent geometry
    mutable GeometryImpl geo_;
    mutable bool builtgeometry_;  //!< true if geometry has been constructed
  }; // end of AlbertaGridEntity codim = 0

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_ENTITY_HH
