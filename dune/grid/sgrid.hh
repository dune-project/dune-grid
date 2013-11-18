// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SGRID_HH
#define DUNE_SGRID_HH

#include <limits>
#include <vector>
#include <stack>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/sgrid/numbering.hh>
#include <dune/grid/common/indexidset.hh>

/*! \file sgrid.hh
   This file documents the DUNE grid interface. We use the special implementation for
   simple structured grid to illustrate the different classes and their members.
 */

namespace Dune {

  //************************************************************************
  /*! define name for floating point type used for coordinates in sgrid.
     You can change the type for coordinates by changing this single typedef.
   */
  typedef double sgrid_ctype;

  // globally define the persistent index type
  const int sgrid_dim_bits = 24;   // bits for encoding each dimension
  const int sgrid_level_bits = 6;  // bits for encoding level number
  const int sgrid_codim_bits = 4;  // bits for encoding codimension

  //************************************************************************
  // forward declaration of templates

  template<int dim, int dimworld, class GridImp> class SGeometry;
  template<int codim, int dim, class GridImp> class SEntity;
  template<int codim, class GridImp> class SEntityPointer;
  template<int codim, class GridImp> class SEntitySeed;
  template<int codim, PartitionIteratorType, class GridImp> class SLevelIterator;
  template<int dim, int dimworld, class ctype> class SGrid;
  template<class GridImp> class SIntersection;
  template<class GridImp> class SIntersectionIterator;
  template<class GridImp> class SHierarchicIterator;

  namespace FacadeOptions
  {

    template<int dim, int dimworld, class ctype, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim,
        SGrid<dim,dimworld,ctype>, SGeometry>
    {
      static const bool v = false;
    };

    template<int dim, int dimworld, class ctype, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim,
        const SGrid<dim,dimworld,ctype>, SGeometry>
    {
      static const bool v = false;
    };

  }

  //************************************************************************
  /*!
     SGeometry realizes the concept of the geometric part of a mesh entity.

     The geometric part of a mesh entity is a \f$d\f$-dimensional object in \f$\mathbf{R}^w\f$
     where \f$d\f$ corresponds the template parameter dim and \f$w\f$ corresponds to the
     template parameter dimworld.

     The \f$d\f$-dimensional object is a polyhedron given by a certain number of corners, which
     are vectors in \f$\mathbf{R}^w\f$.

     The member function global provides a map from a topologically equivalent polyhedron ("reference element")
     in \f$\mathbf{R}^d\f$ to the given polyhedron. This map can be inverted by the member function local, where
     an appropriate projection is applied first, when \f$d\neq w\f$.

     In the case of a structured mesh discretizing a generalized cube this map is linear
     and can be described as \f[ g(l) = s + \sum\limits_{i=0}^{d-1} l_ir^i\f] where \f$s\in\mathbf{R}^w\f$
     is a given position vector, the \f$r^i\in\mathbf{R}^w\f$ are given direction vectors and \f$l\in\mathbf{R}^d\f$
     is a local coordinate within the reference polyhedron. The direction vectors are assumed to be orthogonal
     with respect to the standard Eucliden inner product.

     The \f$d\f$-dimensional reference polyhedron is given
     by the points \f$\{ (x_0,\ldots,x_{d-1}) \ | \ x_i\in\{0,1\}\ \}\f$.

     In order to invert the map for a point \f$p\f$, we have to find a local coordinate \f$l\f$ such
     that \f$g(l)=p\f$. Of course this is only possible if \f$d=w\f$. In the general case \f$d\leq w\f$
     we determine \f$l\f$ such that
     \f[(s,r^k) + \sum\limits_{i=0}^{d-1} l_i (r^i,r^k) = (p,r^k) \ \ \ \forall k=0,\ldots,d-1. \f]

     The resulting system is diagonal since the direction vectors are required to be orthogonal.
   */
  template<int mydim, int cdim, class GridImp>
  class SGeometry
    : public AxisAlignedCubeGeometry<typename GridImp::ctype,mydim,cdim>
  {
  public:
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    /** \brief Set up the geometry
     *
     * \param lower The lower left corner
     * \param A The direction vectors
     *
     * Allows a consistent treatment of all dimensions, including 0 (the vertex).
     */
    void make (const FieldVector<ctype,cdim>& lower,
               const FieldMatrix<ctype,mydim,cdim>& A)
    {
      if (mydim==0) {
        // set up base class
        static_cast< AxisAlignedCubeGeometry<ctype,mydim,cdim> & >( *this ) = AxisAlignedCubeGeometry<ctype,mydim,cdim>(lower);
        return;
      }

      // construct the upper right corner of the cube geometry
      FieldVector<ctype, cdim> upper = lower;
      for (int i=0; i<mydim; i++)
        upper += A[i];

      // look for the directions where the cube is actually extended
      std::bitset<cdim> axes;

      for (size_t i=0; i<cdim; i++)
        if ((upper[i] - lower[i]) > 1e-10)
          axes[i] = true;

      // set up base class
      static_cast< AxisAlignedCubeGeometry<ctype,mydim,cdim> & >( *this ) = AxisAlignedCubeGeometry<ctype,mydim,cdim>(lower, upper, axes);
    }

    //! constructor
    SGeometry ()
      : AxisAlignedCubeGeometry<ctype,mydim,cdim>(FieldVector<ctype,cdim>(0),FieldVector<ctype,cdim>(0))    // anything
    {}
  };


  //************************************************************************
  /*! SEntityBase contains the part of SEntity that can be defined
     without specialization. This is the base for all SEntity classes with dim>0.
   */

  template<int codim, int dim, class GridImp, template<int,int,class> class EntityImp>
  class SEntityBase :
    public EntityDefaultImplementation<codim,dim,GridImp,EntityImp>
  {
    friend class SEntityPointer<codim,GridImp>;
    friend class SIntersectionIterator<GridImp>;
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim< codim >::GeometryImpl GeometryImpl;

  public:
    typedef typename GridImp::ctype ctype;
    typedef typename GridImp::template Codim<codim>::Geometry Geometry;
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! level of this element
    int level () const
    {
      return l;
    }

    //! global index is calculated from the index and grid size
    int globalIndex() const;

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    SEntitySeed<codim, GridImp> seed () const {
      return SEntitySeed<codim, GridImp>(l, index);
    }

    //! return the element type identifier
    GeometryType type () const
    {
      static const GeometryType cubeType(GeometryType::cube,dim-codim);
      return cubeType;
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      if (!builtgeometry) makegeometry();

      // return result
      return Geometry( geo );
    }

    PartitionType partitionType () const { return InteriorEntity; }

    //! constructor
    SEntityBase (GridImp* _grid, int _l, int _index) :
      grid(_grid),
      l(_l),
      index(_index),
      z(grid->z(l,index,codim)),
      builtgeometry(false) {}

    //! empty constructor
    SEntityBase () :
      builtgeometry(false) // mark geometry as not built
    {}

    //! copy constructor
    SEntityBase ( const SEntityBase& other ) :
      grid(other.grid),
      l(other.l),
      index(other.index),
      z(other.z),
      geo(), // do not copy geometry
      builtgeometry(false) // mark geometry as not built
    {}

    //! Reinitialization
    void make (GridImp* _grid, int _l, int _id);

    //! Reinitialization
    void make (int _l, int _id);

    //! geometry of this entity
    void makegeometry () const;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      return grid->persistentIndex(l, codim, z);
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const
    {
      return index;
    }

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      // codim != dim -> there are no copies of entities
      // maxlevel -> ids are fine
      if (codim<dim || l==grid->maxLevel())
        return compressedIndex();

      // this is a vertex which is not on the finest level
      // move coordinates up to maxlevel (multiply by 2 for each level
      array<int,dim> coord;
      for (int k=0; k<dim; k++)
        coord[k] = z[k]*(1<<(grid->maxLevel()-l));

      // compute number with respect to maxLevel
      return grid->n(grid->maxLevel(),coord);
    }

    //! subentity compressed index (not available here)
    int subCompressedIndex (int cd, int i) const
    {
      DUNE_THROW(NotImplemented,"subIndex for entities with codimension > 0 is not implemented");
      return -1;
    }

    //! subentity compressed leaf index (not available here)
    int subCompressedLeafIndex (int cd, int i) const
    {
      DUNE_THROW(NotImplemented,"subIndex for entities with codimension > 0 is not implemented");
      return -1;
    }

  protected:
    // this is how we implement our elements
    GridImp* grid;       //!< grid containes mapper, geometry, etc.
    int l;               //!< level where element is on
    int index;           //!< my consecutive index
    array<int,dim> z;    //!< my coordinate, number of even components = codim
    mutable GeometryImpl geo; //!< geometry, is only built on demand
    mutable bool builtgeometry; //!< true if geometry has been constructed
  };


  /**
     A Grid is a container of grid entities. An entity is parametrized by
     the codimension.  An entity of codimension c in dimension d is a d-c
     dimensional object.
   */
  template<int codim, int dim, class GridImp>
  class SEntity : public SEntityBase<codim,dim,GridImp,SEntity>
  {
    typedef Dune::SEntityBase<codim,dim,GridImp,Dune::SEntity> SEntityBase;
    friend class SEntityPointer<codim,GridImp>;
    friend class SIntersectionIterator<GridImp>;
  public:
    //! constructor
    SEntity (GridImp* _grid, int _l, int _id) :
      SEntityBase(_grid,_l,_id) {}
  };

  /**
     A Grid is a container of grid entities. An entity is parametrized
     by the codimension.  An entity of codimension c in dimension d is a
     d-c dimensional object.

     Entities of codimension 0 ("elements") are defined through template
     specialization. Note that this specialization has an extended
     interface compared to the general case

     Entities of codimension 0 allow to visit all neighbors, where a
     neighbor is an entity of codimension 0 which has a common entity of
     codimension 1 with the entity.  These neighbors are accessed via an
     iterator. This allows the implementation of non-matching meshes. The
     number of neigbors may be different from the number of faces/edges
     of an element!
   */

  /**
     A Grid is a container of grid entities. An entity is parametrized by
     the codimension.  An entity of codimension c in dimension d is a d-c
     dimensional object.

     Entities of codimension=0 ("Cells") are defined through template
     specialization. Note that this specialization has an extended
     interface compared to the general case
   */
  template<int dim, class GridImp>
  class SEntity<0,dim,GridImp> : public SEntityBase<0,dim,GridImp,SEntity>
  {
    enum { dimworld = GridImp::dimensionworld };
    typedef Dune::SEntityBase<0,dim,GridImp,Dune::SEntity> SEntityBase;
    using SEntityBase::grid;
    using SEntityBase::l;
    using SEntityBase::index;
    using SEntityBase::z;

    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 0 >::LocalGeometryImpl LocalGeometryImpl;

    friend class SEntityPointer<0,GridImp>;
    friend class SIntersectionIterator<GridImp>;

  public:
    typedef typename GridImp::ctype ctype;
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;
    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridImp::HierarchicIterator HierarchicIterator;
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! make HierarchicIterator a friend
    friend class SHierarchicIterator<GridImp>;

    /**
       Intra-element access to entities of codimension cc > codim.
       Return number of entities with codimension cc.
     */
    template<int cc> int count () const;

    /**
       Provide access to mesh entity i of given codimension. Entities
       are numbered 0 ... count<cc>()-1
     */
    template<int cc> typename Codim<cc>::EntityPointer subEntity (int i) const;

    //! subentity compressed index
    int subCompressedIndex (int codim, int i) const
    {
      if (codim==0) return this->compressedIndex();
      // compute subIndex
      return (this->grid)->n(this->l, this->grid->subz(this->z,i,codim));
    }

    /*! subentity leaf index
       \todo add handling of not-leaf vertices
     */
    int subCompressedLeafIndex (int codim, int i) const
    {
      if (codim==0) return this->compressedLeafIndex();

      assert(this->l == this->grid->maxLevel());
      // compute subIndex
      return (this->grid)->n(this->l, this->grid->subz(this->z,i,codim));
    }

    //! subentity persistent index
    PersistentIndexType subPersistentIndex (int codim, int i) const
    {
      if (codim==0) return this->persistentIndex();
      // compute subId
      return this->grid->persistentIndex(this->l, codim, this->grid->subz(this->z,i,codim));
    }

    /**
       Intra-level access to intersections with neighboring elements.  A
       neighbor is an entity of codimension 0 which has an entity of
       codimension 1 in commen with this entity. Access to neighbors is
       provided using iterators. This allows meshes to be
       nonmatching. Returns iterator referencing the first neighbor.
     */
    IntersectionIterator ibegin () const;
    IntersectionIterator ileafbegin () const;
    IntersectionIterator ilevelbegin () const;
    //! Reference to one past the last intersection
    IntersectionIterator iend () const;
    IntersectionIterator ileafend () const;
    IntersectionIterator ilevelend () const;

    /**
       @brief Inter-level access to father element on coarser grid.

       Assumes that meshes are nested.
     */
    EntityPointer father () const;

    //! returns true if father entity exists
    bool hasFather () const
    {
      return (this->level()>0);
    }

    //! return true if the entity is leaf
    bool isLeaf () const
    {
      return ( this->grid->maxLevel() == this->level() );
    }

    /**
       @brief Location of this element relative to the reference element element of the father.

       This is sufficient to interpolate all dofs in conforming case.
       Nonconforming may require access to neighbors of father and
       computations with local coordinates.  On the fly case is somewhat
       inefficient since dofs are visited several times.  If we store
       interpolation matrices, this is tolerable. We assume that
       on-the-fly implementation of numerical algorithms is only done for
       simple discretizations.  Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather () const;

    /**
       @brief Inter-level access to son elements on higher levels<=maxLevel.

       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    HierarchicIterator hbegin (int maxLevel) const;

    //! Returns iterator to one past the last son
    HierarchicIterator hend (int maxLevel) const;

    // members specific to SEntity
    //! constructor
    SEntity (GridImp* _grid, int _l, int _index) :
      SEntityBase(_grid,_l,_index),
      built_father(false)
    {}

    SEntity (const SEntity& other ) :
      SEntityBase(other.grid, other.l, other.index ),
      built_father(false)
    {}

    //! Reinitialization
    void make (GridImp* _grid, int _l, int _id)
    {
      SEntityBase::make(_grid,_l,_id);
      built_father = false;
    }

    //! Reinitialization
    void make (int _l, int _id)
    {
      SEntityBase::make(_l,_id);
      built_father = false;
    }

  private:

    SEntity();

    mutable bool built_father;
    mutable int father_index;
    mutable LocalGeometryImpl in_father_local;
    void make_father() const;
  };


  //************************************************************************
  /*! Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HIerarchicIterator,
     starting from a given entity.
     This is redundant but important for memory efficient implementations of unstructured
     hierarchically refined meshes.
   */
  struct SHierarchicStackElem {
    int l;
    int index;
    SHierarchicStackElem () : l(-1), index(-1) {}
    SHierarchicStackElem (int _l, int _index) {l=_l; index=_index;}
    bool operator== (const SHierarchicStackElem& s) const {return !operator!=(s);}
    bool operator!= (const SHierarchicStackElem& s) const {return l!=s.l || index!=s.index;}
  };

  template<class GridImp>
  class SHierarchicIterator :
    public Dune::SEntityPointer <0,GridImp>
  {
    friend class SHierarchicIterator<const GridImp>;
    enum { dim = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };
    typedef Dune::SEntityPointer<0,GridImp> SEntityPointer;
    using SEntityPointer::realEntity;
    using SEntityPointer::grid;
    using SEntityPointer::l;
    using SEntityPointer::index;
  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ctype;

    //! increment
    void increment();

    /*! constructor. Here is how it works: If with_sons is true, push start
       element and all its sons on the stack, so the initial element is popped
       last. For an end iterator, push the starting element and no sons. Then
       the iteration will stop when both iterators have the same id AND the
       stack is empty
     */
    SHierarchicIterator (GridImp* _grid,
                         const Dune::SEntity<0,GridImp::dimension,GridImp>& _e,
                         int _maxLevel, bool makeend) :
      SEntityPointer(_grid,_e.level(),_e.compressedIndex())
    {
      // without sons, we are done
      // (the end iterator is equal to the calling iterator)
      if (makeend) return;

      // remember element where begin has been called
      orig_l = this->entity().level();
      orig_index = _grid->getRealImplementation(this->entity()).compressedIndex();

      // push original element on stack
      SHierarchicStackElem originalElement(orig_l, orig_index);
      stack.push(originalElement);

      // compute maxLevel
      maxLevel = std::min(_maxLevel,this->grid->maxLevel());

      // ok, push all the sons as well
      push_sons(orig_l,orig_index);

      // and pop the first son
      increment();
    }

  private:
    int maxLevel;              //!< maximum level of elements to be processed
    int orig_l, orig_index;       //!< element where begin was called (the root of the tree to be processed)

    //!< stack holding elements to be processed
    std::stack<SHierarchicStackElem, Dune::ReservedVector<SHierarchicStackElem,GridImp::MAXL> > stack;

    void push_sons (int level, int fatherid); //!< push all sons of this element on the stack
  };

  //************************************************************************
  /*! Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1 with the entity.
     These neighbors are accessed via a IntersectionIterator. This allows the implementation of
     non-matching meshes. The number of neigbors may be different from the number of faces/edges
     of an element!
   */
  template<class GridImp>
  class SIntersectionIterator
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

    friend class SIntersection<GridImp>;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef Dune::SIntersection<GridImp> IntersectionImp;
    typedef Dune::Intersection< const GridImp, Dune::SIntersection< const GridImp > > Intersection;
    //! know your own dimension
    enum { dimension=dim };
    //! know your own dimension of world
    enum { dimensionworld=dimworld };
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! equality
    bool equals(const SIntersectionIterator<GridImp>& i) const;
    //! increment
    void increment();

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return intersection;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const;

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const;

    //! return true if intersection is with boundary.
    bool boundary () const;

    int boundaryId () const {
      if (boundary()) return count + 1;
      return 0;
    }

    int boundarySegmentIndex () const {
      if (boundary())
        return grid->boundarySegmentIndex(self.level(), count, zred);
      return -1;
    }

    //! return true if neighbor on this level exists
    bool neighbor () const;

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of the element where iteration started.
     */
    LocalGeometry geometryInInside () const;
    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of neighbor
     */
    LocalGeometry geometryInOutside () const;
    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in GLOBAL coordinates of the element where iteration started.
     */
    Geometry geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      static const GeometryType cubeType(GeometryType::cube,dim-1);
      return cubeType;
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const;
    //! local index of codim 1 entity in neighbor where intersection is contained in
    int indexInOutside () const;

    //! constructor
    SIntersectionIterator (GridImp* _grid, const SEntity<0,dim,GridImp >* _self, int _count) :
      self(*_self), ne(self), grid(_grid),
      partition(_grid->partition(grid->getRealImplementation(ne).l,_self->z)),
      zred(_grid->compress(grid->getRealImplementation(ne).l,_self->z)),
      intersection(IntersectionImp(*this))
    {
      // make neighbor
      make(_count);
    }

    SIntersectionIterator (const SIntersectionIterator & other) :
      self(other.self), ne(other.ne), grid(other.grid),
      partition(other.partition), zred(other.zred),
      count(other.count), valid_count(other.valid_count),
      valid_nb(other.valid_nb), is_on_boundary(other.is_on_boundary),
      built_intersections(false),
      intersection(IntersectionImp(*this))
    {}

    //! assignment operator
    SIntersectionIterator& operator = (const SIntersectionIterator& other)
    {
      /* We can't assign the grid */
      assert(grid == other.grid);

      /* Assign data from other */
      self = other.self;
      ne = other.ne;
      partition = other.partition;
      zred = other.zred;
      count = other.count;
      valid_count = other.valid_count;
      valid_nb = other.valid_nb;
      is_on_boundary = other.is_on_boundary;

      /* mark cached data as invalid */
      built_intersections = false;

      return *this;
    }

  private:
    void make (int _count) const;         //!< reinitialze iterator with given neighbor
    void makeintersections () const;      //!< compute intersections
    EntityPointer self;                   //!< EntityPointer for myself
    mutable EntityPointer ne;             //!< EntityPointer for neighbor
    const GridImp * grid;                 //!< Pointer to the grid
    int partition;                        //!< partition number of self, needed for coordinate expansion
    array<int,dim> zred;                  //!< reduced coordinates of myself, allows easy computation of neighbors
    mutable int count;                    //!< number of neighbor
    mutable bool valid_count;             //!< true if count is in range
    mutable bool valid_nb;                //!< true if nb is initialized
    mutable bool is_on_boundary;          //!< true if neighbor is otside the domain
    mutable bool built_intersections;     //!< true if all intersections have been built
    mutable LocalGeometryImpl is_self_local;  //!< intersection in own local coordinates
    mutable GeometryImpl is_global;           //!< intersection in global coordinates, map consistent with is_self_local
    mutable LocalGeometryImpl is_nb_local;    //!< intersection in neighbors local coordinates
    Intersection intersection;
  };

  template<class GridImp>
  class SIntersection
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename Geometry::LocalCoordinate LocalCoordinate;
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef Dune::Intersection< const GridImp, Dune::SIntersectionIterator< const GridImp > > Intersection;
    //! know your own dimension
    enum { dimension=dim };
    //! know your own dimension of world
    enum { dimensionworld=dimworld };
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    bool boundary () const
    {
      return is.boundary();
    }

    /*! @brief Identifier for boundary segment from macro grid. */
    int boundaryId () const
    {
      return is.boundaryId();
    }

    /*! @brief index of the boundary segment within the macro grid  */
    size_t boundarySegmentIndex () const
    {
      return is.boundarySegmentIndex();
    }

    /*! @brief return true if intersection is shared with another element. */
    bool neighbor () const
    {
      return is.neighbor();
    }

    /*! @brief return EntityPointer to the Entity on the inside of this intersection. */
    EntityPointer inside() const
    {
      return is.inside();
    }

    /*! @brief return EntityPointer to the Entity on the outside of this intersection. */
    EntityPointer outside() const
    {
      return is.outside();
    }

    /*! @brief return true if intersection is conform. */
    bool conforming () const
    {
      return true;
    }

    /*! @brief geometrical information about this intersection in local coordinates of the inside() entity. */
    LocalGeometry geometryInInside () const
    {
      return is.geometryInInside();
    }

    /*! @brief geometrical information about this intersection in local coordinates of the outside() entity. */
    LocalGeometry geometryInOutside () const
    {
      return is.geometryInOutside();
    }

    /*! @brief geometrical information about the intersection in global coordinates. */
    Geometry geometry () const
    {
      return is.geometry();
    }

    /*! @brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return is.type();
    }

    /*! @brief Local index of codim 1 entity in the inside() entity where intersection is contained in */
    int indexInInside () const
    {
      return is.indexInInside();
    }

    /*! @brief Local index of codim 1 entity in outside() entity where intersection is contained in */
    int indexInOutside () const
    {
      return is.indexInOutside();
    }

    /*! @brief Return an outer normal (length not necessarily 1) */
    GlobalCoordinate outerNormal (const LocalCoordinate& local) const
    {
      return centerUnitOuterNormal();
    }

    /*! @brief return outer normal scaled with the integration element */
    GlobalCoordinate integrationOuterNormal (const LocalCoordinate& local) const
    {
      FieldVector<ctype, dimworld> n = centerUnitOuterNormal();
      n *= is.geometry().integrationElement(local);
      return n;
    }

    /*! @brief Return unit outer normal (length == 1)  */
    GlobalCoordinate unitOuterNormal (const LocalCoordinate& local) const
    {
      return centerUnitOuterNormal();
    }

    /*! @brief Return unit outer normal (length == 1) */
    GlobalCoordinate centerUnitOuterNormal () const
    {
      FieldVector<ctype, dimworld> normal(0.0);
      normal[is.count/2] =  (is.count%2) ? 1.0 : -1.0;
      return normal;
    }

    //! constructor
    SIntersection (const SIntersectionIterator<GridImp> & is_) : is(is_) {}

  private:
#ifndef DOXYGEN // doxygen can't handle this recursive usage
    const SIntersectionIterator<GridImp> & is;
#endif
  };

  //************************************************************************

  /** \brief a stack of pointers with auto destruction if the stack is
      destructed
   */
  template <class T>
  class AutoPtrStack : public std::stack<T*>
  {
  public:
    ~AutoPtrStack()
    {
      while(! this->empty())
      {
        T* e = this->top();
        delete e;
        this->pop();
      }
    }
  };

  /*! Acts as a pointer to an  entities of a given codimension.
   */
  template<int codim, class GridImp>
  class SEntityPointer
  {
    enum { dim = GridImp::dimension };
    friend class SIntersectionIterator<GridImp>;
  public:
    typedef SEntityPointer<codim,GridImp> EntityPointerImp;
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    //! codimension of entity pointer
    enum { codimension = codim };

    //! equality
    bool equals(const SEntityPointer<codim,GridImp>& i) const;
    //! dereferencing
    Entity& dereference() const;
    //! ask for level of entity
    int level () const;

    //! constructor
    SEntityPointer (GridImp * _grid, int _l, int _index) :
      grid(_grid), l(_l), index(_index),
      e(0)
    {}

    //! constructor
    SEntityPointer (const SEntity<codim,dim,GridImp> & _e) :
      grid(_e.grid), l(_e.l), index(_e.index),
      e(0)
    {}

    //! constructor
    SEntityPointer (const SEntityPointer<codim,GridImp>& other) :
      grid(other.grid), l(other.l), index(other.index),
      e( 0 )
    {}

    //! destructor pointer
    ~SEntityPointer()
    {
      if( e )
        enStack().push( e );
#ifndef NDEBUG
      index = -1;
#endif
    }

    //! assignment operator
    SEntityPointer& operator = (const SEntityPointer& other)
    {
      grid = other.grid;
      l = other.l;
      index = other.index;

      // free current entity
      if( e )
        enStack().push( e );
      e = 0;

      return *this;
    }

  protected:
    inline SEntity<codim,dim,GridImp>& realEntity() const
    {
      return grid->getRealImplementation(entity());
    }

    inline Entity& entity() const
    {
      if( ! e )
      {
        e = getEntity( grid, l, index );
      }
      return *e;
    }

    typedef AutoPtrStack< Entity > EntityStackType;
    static inline EntityStackType& enStack()
    {
      static EntityStackType eStack;
      return eStack;
    }

    inline Entity* getEntity(GridImp* _grid, int _l, int _id ) const
    {
      // get stack reference
      EntityStackType& enSt = enStack();

      if( enSt.empty() )
      {
        return (new Entity(SEntity<codim,dim,GridImp>(_grid, _l, _id)));
      }
      else
      {
        Entity* e = enSt.top();
        enSt.pop();
        grid->getRealImplementation(*e).make(_grid, _l,_id);
        return e;
      }
    }

    GridImp* grid;               //!< my grid
    int l;                       //!< level where element is on
    mutable int index;           //!< my consecutive index
    mutable Entity* e;           //!< virtual entity
  };

  /*! describes the minimal information necessary to create a fully functional SEntity
   */
  template<int codim, class GridImp>
  class SEntitySeed
  {
    enum { dim = GridImp::dimension };
  public:
    enum { codimension = codim };

    //! default constructor (invalid)
    SEntitySeed () :
      _l(-1), _index(0)
    {}

    //! constructor
    SEntitySeed (int l, int index) :
      _l(l), _index(index)
    {}

    //! check whether the EntitySeed refers to a valid Entity
    bool isValid() const
    {
      return _l != -1;
    }

    int level () const { return this->_l; }
    int index () const { return this->_index; }

  private:
    int _l;                       //!< level where element is on
    int _index;                   //!< my consecutive index
  };

  //************************************************************************


  /*! Enables iteration over all entities of a given codimension and level of a grid.
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class SLevelIterator :
    public Dune::SEntityPointer <codim,GridImp>
  {
    friend class SLevelIterator<codim, pitype,const GridImp>;
    enum { dim = GridImp::dimension };
    typedef Dune::SEntityPointer<codim,GridImp> SEntityPointer;
    using SEntityPointer::realEntity;
    using SEntityPointer::l;
    using SEntityPointer::index;
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! increment
    void increment();

    //! constructor
    SLevelIterator (GridImp * _grid, int _l, int _id) :
      SEntityPointer(_grid,_l,_id) {}
  };


  //========================================================================
  /*!
     \brief implementation of index set

   */
  //========================================================================

  template<class GridImp>
  class SGridLevelIndexSet : public IndexSet<GridImp,SGridLevelIndexSet<GridImp> >
  {
    typedef SGridLevelIndexSet< GridImp > This;
    typedef IndexSet< GridImp, This > Base;

    enum { dim = GridImp::dimension };

  public:

    //! constructor stores reference to a grid and level
    SGridLevelIndexSet ( const GridImp &g, int l )
      : grid( g ),
        level( l )
    {
      // TODO move list of geometrytypes to grid, can be computed static (singleton)
      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryType(GeometryType::cube,GridImp::dimension-codim));
    }

    //! get index of an entity
    template<int cd>
    int index (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return grid.getRealImplementation(e).compressedIndex();
    }

    template< int cc >
    int subIndex ( const typename GridImp::Traits::template Codim< cc >::Entity &e,
                   int i, unsigned int codim ) const
    {
      if( cc == 0 )
        return grid.getRealImplementation(e).subCompressedIndex(codim, i);
      else
        DUNE_THROW( NotImplemented, "subIndex for higher codimension entity not implemented for SGrid." );
    }

    // return true if the given entity is contained in \f$E\f$.
    template< class EntityType >
    bool contains ( const EntityType &e ) const
    {
      return (e.level() == level);
    }

    //! get number of entities of given type and level (the level is known to the object)
    int size (GeometryType type) const
    {
      return grid.size( level, type );
    }

    //! return size of set for a given codim
    int size (int codim) const
    {
      return grid.size( level, codim );
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return mytypes[codim];
    }

  private:
    const GridImp& grid;
    int level;
    std::vector<GeometryType> mytypes[GridImp::dimension+1];
  };



  //========================================================================
  /*!
     \brief persistent, globally unique Ids

   */
  //========================================================================

  template<class GridImp>
  class SGridGlobalIdSet :
    public IdSet<GridImp,SGridGlobalIdSet<GridImp>, typename remove_const<GridImp>::type::PersistentIndexType>
    /*
       We used the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
  {
    typedef SGridGlobalIdSet< GridImp > This;

  public:

    //! define the type used for persistent indices
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    typedef typename remove_const<GridImp>::type::PersistentIndexType IdType;

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return GridImp::getRealImplementation(e).persistentIndex();
    }

    //! get id of subentity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    IdType subId ( const typename remove_const< GridImp >::type::Traits::template Codim< 0 >::Entity &e,
                   int i, unsigned int codim ) const
    {
      return GridImp::getRealImplementation(e).subPersistentIndex(codim, i);
    }
  };


  template<int dim, int dimworld, class ctype>
  struct SGridFamily
  {
    typedef GridTraits<dim,dimworld,Dune::SGrid<dim,dimworld,ctype>,
        SGeometry,SEntity,
        SEntityPointer,SLevelIterator,
        SIntersection,              // leaf intersection
        SIntersection,              // level intersection
        SIntersectionIterator,              // leaf  intersection iter
        SIntersectionIterator,              // level intersection iter
        SHierarchicIterator,
        SLevelIterator,
        SGridLevelIndexSet<const SGrid<dim,dimworld,ctype> >,
        SGridLevelIndexSet<const SGrid<dim,dimworld,ctype> >,
        SGridGlobalIdSet<const SGrid<dim,dimworld,ctype> >,
        bigunsignedint<dim*sgrid_dim_bits+sgrid_level_bits+sgrid_codim_bits>,
        SGridGlobalIdSet<const SGrid<dim,dimworld,ctype> >,
        bigunsignedint<dim*sgrid_dim_bits+sgrid_level_bits+sgrid_codim_bits>,
        CollectiveCommunication<Dune::SGrid<dim,dimworld,ctype> >,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        SEntitySeed>
    Traits;
  };


  //************************************************************************
  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief A structured mesh in d dimensions consisting of "cubes" (pilot implementation of the %Dune grid interface, for debugging only).
     \ingroup GridImplementations
     \ingroup SGrid

          This module describes the pilot implementation of the %Dune grid interface.
          It implements the grid interface for simple structured meshes.

          \warning SGrid is slow. It is intended for debugging only.

          The following class diagram shows how the classes are related with
          each other:

          \image html sgridclasses.png "Class diagram for classes in the grid interface"
          \image latex sgridclasses.eps "Class diagram for classes in the grid interface" width=\textwidth

          Short description of the classes:

          - SGeometry is a class template providing the geometric part of a grid entity, i.e. a general polyhedron
          with a mapping from a reference polyhedron to the actual polyhedron.

          - SLevelIterator is a class template which allows to iterate over all grid entities of a given
          codimension and level.

          - SEntity is a class template realizing the grid entities. Grid entities are the constituents
          of a grid. Grid entities of codimension 0 and codimension dim are defines through specialization.
          Entities can be used as template parameters to generic algorithms. Each entity must therefore
          provide the nested classes Geometry, LevelIterator, HierarchicIterator and IntersectionIterator.
          Geometry and LevelIterator are derived from the classes SELement and SLevelIterator.
          Note that entities of codimension 0 and dim have an extended interface.

          - SEntity::IntersectionIterator provides access to all entities of codimension 0 sharing an object of codimension 1
          with the given entity of codimension 0. This interface covers nonmatching grids.

          - SEntity::HierarchicIterator provides access to the sons of an entity of codimension 0.

          - SGrid is conceptualized as a container of grid entities of various codimensions. Since grids
          are used as template parameters to generic algorithms they must include the nested classes
          LevelIterator and Entity which are derived from SLevelIterator and SEntity.

     A Grid is a container of grid entities. Given a dimension dim these entities have a
     codimension codim with 0 <= codim <= dim.

     The Grid is assumed to be hierachically refined and nested. It enables iteration over
     entities of a given level and codimension.

     All information is provided to allocate degrees of freedom in appropriate vector
     data structures.

     \note When SGrid is instantiated with dimworld strictly greater than dim, the result is a
     dim-dimensional structured grid which is embedded in the first dim components of
     dimworld-dimensional Euclidean space.
   */
  template<int dim, int dimworld, typename _ctype = sgrid_ctype>
  class SGrid : public GridDefaultImplementation <dim,dimworld,_ctype,SGridFamily<dim,dimworld,_ctype> >
  {
  public:
    typedef SGridFamily<dim,dimworld,_ctype> GridFamily;
    typedef bigunsignedint<dim*sgrid_dim_bits+sgrid_level_bits+sgrid_codim_bits> PersistentIndexType;

    // need for friend declarations in entity
    typedef SGridLevelIndexSet<SGrid<dim,dimworld> > LevelIndexSetType;
    typedef SGridLevelIndexSet<SGrid<dim,dimworld> > LeafIndexSetType;
    typedef SGridGlobalIdSet<SGrid<dim,dimworld> > GlobalIdSetType;

    typedef typename SGridFamily<dim,dimworld,_ctype>::Traits Traits;

    //! maximum number of levels allowed
    enum { MAXL=32 };

    //! define type used for coordinates in grid module
    typedef _ctype ctype;

    // constructors

    /*! @brief Make an SGrid from extend and number of cells per direction

       \param[in] N_ number of cells in each direction on coarsest level
       \param[in] H_ extend of the unit cube in each dimension

       Note: The origin of the cube is always at (0,0,...,0), only the extend is given.
     */
    SGrid (const int * const N_, const ctype * const H_);

    /*! @brief Make an SGrid from position, extend and number of cells per direction

       \param[in] N_ number of cells in each direction on coarsest level
       \param[in] L_ position of origin of the cube
       \param[in] H_ position of the upper right corner of the cube

     */
    SGrid (const int * const N_, const ctype * const L_, const ctype * const H_);

    /*! @brief Make an SGrid from position, extend and number of cells per direction

       \param[in] N_ number of cells in each direction on coarsest level
       \param[in] L_ position of origin of the cube
       \param[in] H_ position of the upper right corner of the cube

       Note: This constructor uses FieldVectors instead of built-in arrays. This is compatible
          with the YaspGrid class.
     */
    SGrid (FieldVector<int,dim> N_, FieldVector<ctype,dim> L_, FieldVector<ctype,dim> H_);

    //! @brief empty constructor making grid of unit square discretized with one cell
    SGrid ();

    //! @brief SGrid destructor
    ~SGrid ();

    /*! Return maximum level defined in this grid. Levels are numbered
          0 ... maxLevel with 0 the coarsest level.   */
    int maxLevel() const;

    //! Iterator to first entity of given codim on level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const;

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      return lbegin<cd,All_Partition>(level);
    }

    //! one past the end on this level
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      return lend<cd,All_Partition>(level);
    }

    //! return LeafIterator which points to the first entity
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const;

    //! one past the end on the leaf level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const;

    //! return LeafIterator which points to the first entity
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      return leafbegin<cd,All_Partition>();
    }

    //! return LeafIterator which points behind the last entity
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      return leafend<cd,All_Partition>();
    }

    // \brief obtain EntityPointer from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::EntityPointer
    entityPointer(const Seed& seed) const
    {
      enum { codim = Seed::codimension };
      return SEntityPointer<codim,const SGrid<dim,dimworld> >(this,
                                                              this->getRealImplementation(seed).level(),
                                                              this->getRealImplementation(seed).index());
    }

    /*! The communication interface
          @tparam T array class holding data associated with the entities
          @tparam P type used to gather/scatter data in and out of the message buffer
          @tparam codim communicate entites of given codim
          @param t array holding data associated with the entities
          @param iftype one of the predefined interface types, throws error if it is not implemented
          @param dir choose beetween forward and backward communication
          @param level communicate for entities on the given level

          Implements a generic communication function sending an object of type P for each entity
       in the intersection of two processors. P has two methods gather and scatter that implement
       the protocol. Therefore P is called the "protocol class".
     */
    template<class T, template<class> class P, int codim>
    void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      // SGrid is sequential and has no periodic boundaries, so do nothing ...
      return;
    }

    //! number of grid entities per level and codim
    int size (int level, int codim) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return size(maxLevel(),codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? size(level,dim-type.dim()) : 0;
    }

    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return size(maxLevel(),type);
    }

    //! \brief returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const
    {
      return boundarysize;
    }

    //! number of grid entities of all level for given codim
    int global_size (int codim) const;

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim)
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim)
    {
      return 0;
    }

    // these are all members specific to sgrid

    /** \brief Refine mesh globally by one refCount levels */
    void globalRefine (int refCount);

    /** \brief Get number of elements in each coordinate direction */
    const array<int, dim>& dims(int level) const {
      return N[level];
    }

    /** \brief Get lower left corner */
    const FieldVector<ctype, dimworld>& lowerLeft() const {
      return low;
    }

    /** \brief Get upper right corner */
    FieldVector<ctype, dimworld> upperRight() const {
      return H;
    }

    //! map adapt to global refine
    bool adapt ()
    {
      globalRefine(1);
      return true;
    }

    // The new index sets from DDM 11.07.2005
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LocalIdSet& localIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      assert(level>=0 && level<=maxLevel());
      return *(indexsets[level]);
    }

    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return *indexsets.back();
    }

    /*!
       @name dummy parallel functions
       @{
     */
    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {}

    const CollectiveCommunication<SGrid>& comm () const
    {
      return ccobj;
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int codim) const
    {
      return 0;
    }

    /*
       @}
     */

  private:
    /*
       Make associated classes friends to grant access to the real entity
     */
    friend class Dune::SGridLevelIndexSet<Dune::SGrid<dim,dimworld> >;
    friend class Dune::SGridGlobalIdSet<Dune::SGrid<dim,dimworld> >;
    friend class Dune::SIntersectionIterator<Dune::SGrid<dim,dimworld> >;
    friend class Dune::SHierarchicIterator<Dune::SGrid<dim,dimworld> >;
    friend class Dune::SEntity<0,dim,Dune::SGrid<dim,dimworld> >;

    friend class Dune::SGridLevelIndexSet<const Dune::SGrid<dim,dimworld> >;
    friend class Dune::SGridGlobalIdSet<const Dune::SGrid<dim,dimworld> >;
    friend class Dune::SIntersectionIterator<const Dune::SGrid<dim,dimworld> >;
    friend class Dune::SHierarchicIterator<const Dune::SGrid<dim,dimworld> >;
    friend class Dune::SEntity<0,dim,const Dune::SGrid<dim,dimworld> >;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Dune::SEntityBase;

    template<int codim_, class GridImp_>
    friend class Dune::SEntityPointer;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    //! map expanded coordinates to position
    FieldVector<ctype, dimworld> pos (int level, array<int,dim>& z) const;

    //! compute codim from coordinate
    int calc_codim (int level, const array<int,dim>& z) const;

    //! compute number from expanded coordinate
    int n (int level, const array<int,dim>& z) const;

    //! compute coordinates from number and codimension
    array<int,dim> z (int level, int i, int codim) const;

    //! compute zentity of subentity of given codim
    array<int,dim> subz (const array<int,dim> & z, int i, int codim) const;

    //! compress from expanded coordinates to grid for a single partition number
    array<int,dim> compress (int level, const array<int,dim>& z) const;

    //! expand with respect to partition number
    array<int,dim> expand (int level, const array<int,dim>& r, int b) const;

    /*! There are \f$2^d\f$ possibilities of having even/odd coordinates.
          The binary representation is called partition number.
     */
    int partition (int level, const array<int,dim>& z) const;

    //! given reduced coordinates of an element, determine if element is in the grid
    bool exists (int level, const array<int,dim>& zred) const;

    // compute boundary segment index for a given zentity and a face
    int boundarySegmentIndex (int l, int face, const array<int,dim> & zentity) const
    {
      array<int,dim-1> zface;
      int dir = face/2;
      int side = face%2;
      // compute z inside the global face
      for (int i=0; i<dir; i++) zface[i] = zentity[i]/(1<<l);
      for (int i=dir+1; i<dim; i++) zface[i-1] = zentity[i]/(1<<l);
      zface = boundarymapper[dir].expand(zface, 0);
      // compute index in the face
      int index = boundarymapper[dir].n(zface);
      // compute offset
      for (int i=0; i<dir; i++)
        index += 2*boundarymapper[i].elements(0);
      index += side*boundarymapper[dir].elements(0);
      return index;
    }

    // compute persistent index for a given zentity
    PersistentIndexType persistentIndex (int l, int codim, const array<int,dim> & zentity) const
    {
      if (codim!=dim)
      {
        // encode codim, this would actually not be necessary
        // because z is unique in codim
        PersistentIndexType id(codim);

        // encode level
        id = id << sgrid_level_bits;
        id = id+PersistentIndexType(l);

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << sgrid_dim_bits;
          id = id+PersistentIndexType(zentity[i]);
        }

        return id;
      }
      else
      {
        // determine min number of trailing zeroes
        // consider that z is on the doubled grid !
        int trailing = 1000;
        for (int i=0; i<dim; i++)
        {
          // count trailing zeros
          int zeros = 0;
          for (int j=0; j<l; j++)
            if (zentity[i]&(1<<(j+1)))
              break;
            else
              zeros++;
          trailing = std::min(trailing,zeros);
        }

        // determine the level of this vertex
        int level = l-trailing;

        // encode codim
        PersistentIndexType id(dim);

        // encode level
        id = id << sgrid_level_bits;
        id = id+PersistentIndexType(level);

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << sgrid_dim_bits;
          id = id+PersistentIndexType(zentity[i]>>trailing);
        }

        return id;
      }
    }

    // disable copy and assign
    SGrid(const SGrid &) {}
    SGrid & operator = (const SGrid &) { return *this; }
    // generate SGrid
    void makeSGrid (const array<int,dim>& N_, const FieldVector<ctype, dim>& L_, const FieldVector<ctype, dim>& H_);

    /*
       internal data
     */
    CollectiveCommunication<SGrid> ccobj;

    ReservedVector<SGridLevelIndexSet<const SGrid<dim,dimworld> >*, MAXL> indexsets;
    SGridGlobalIdSet<const SGrid<dim,dimworld> > theglobalidset;

    int L;                        // number of levels in hierarchic mesh 0<=level<L
    FieldVector<ctype, dim> low;  // lower left corner of the grid
    FieldVector<ctype, dim> H;    // length of cube per direction
    std::vector<array<int,dim> > N;            // number of elements per direction for each level
    std::vector<FieldVector<ctype, dim> > h;   // mesh size per direction for each level
    mutable CubeMapper<dim> *mapper; // a mapper for each level

    // boundary segement index set
    array<CubeMapper<dim-1>, dim> boundarymapper; // a mapper for each coarse grid face
    int boundarysize;
  };

  namespace Capabilities
  {

    /** \struct isParallel
       \ingroup SGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup SGrid
     */

    /** \brief SGrid has only one geometry type for codim 0 entities
       \ingroup SGrid
     */
    template<int dim, int dimw>
    struct hasSingleGeometryType< SGrid<dim,dimw> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< dim > :: type :: id ;
    };

    /** \brief SGrid is a Cartesian grid
        \ingroup SGrid
     */
    template<int dim, int dimw>
    struct isCartesian< SGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief SGrid has entities for all codimension
       \ingroup SGrid
     */
    template<int dim, int dimw, int cdim>
    struct hasEntity< SGrid<dim,dimw>, cdim>
    {
      static const bool v = true;
    };

    /** \brief SGrid is levelwise conforming
       \ingroup SGrid
     */
    template<int dim, int dimw>
    struct isLevelwiseConforming< SGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief SGrid is leafwise conforming
       \ingroup SGrid
     */
    template<int dim, int dimw>
    struct isLeafwiseConforming< SGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

} // end namespace Dune

#include "sgrid/sgrid.cc"

#endif
