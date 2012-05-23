// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_HH
#define DUNE_ALU2DGRID_INTERSECTION_HH

#include <stack>
#include <utility>

#include <dune/common/nullptr.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid/alugrid/2d/alu2dinclude.hh>
#include <dune/grid/alugrid/2d/entity.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template<int cd, int dim, class Grid>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class Grid >
  class ALU2dGridLevelIterator;
  template<int cd, class Grid >
  class ALU2dGridEntityPointer;
  template<int mydim, int coorddim, class Grid>
  class ALU2dGridGeometry;
  template<class Grid>
  class ALU2dGridHierarchicIterator;
  template<int codim, PartitionIteratorType pitype, class Grid>
  class ALU2dGridLeafIterator;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< class Grid, class Info >
  class ALU2dGridIntersectionBase;

  template< class Grid >
  class ALU2dGridLeafIntersection;

  template< class Grid >
  class ALU2dGridLevelIntersection;

  template< class Grid >
  class ALU2dGridLeafIntersectionIterator;

  template< class Grid >
  class ALU2dGridLevelIntersectionIterator;



  // ALU2DIntersectionGeometryStorage
  // --------------------------------

  template< class LocalGeometryImpl >
  class ALU2DIntersectionGeometryStorage
  {
    typedef ALU2DIntersectionGeometryStorage< LocalGeometryImpl > This;

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
    static const This &instance ()
    {
      static const This geomStorage;
      return geomStorage;
    }
  };



  // ALU2dGridIntersectionBase
  // -------------------------

  template< class Grid, class Info >
  class ALU2dGridIntersectionBase
  {
    typedef ALU2dGridIntersectionBase< Grid, Info > This;

    friend class ALU2dGridLeafIntersectionIterator< Grid >;
    friend class ALU2dGridLevelIntersectionIterator< Grid >;

    static const ALU2DSPACE ElementType eltype = Grid::elementType;

    typedef typename Grid::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename Grid::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

  public:
    typedef Info IntersectionInfo;

    typedef alu2d_ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld  = Grid::dimensionworld;

    typedef typename Grid::GridObjectFactoryType Factory;

    typedef typename Grid::template Codim<0>::Entity Entity;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

    typedef typename Grid::template Codim<1>::Geometry Geometry;
    typedef typename Grid::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dimension,Grid> EntityImp;

    typedef FieldVector< ctype, dimensionworld > NormalType;
    typedef FieldVector< ctype, dimension-1 > LocalCoordinate;

    typedef ALU2dGridEntityPointer<0,Grid> EntityPointerImp;

    typedef typename ALU2dImplTraits< dimensionworld, eltype >::ThinelementType ThinelementType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::HElementType HElementType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::HBndElType HBndElType;

    // type of local geometry storage
    typedef ALU2DIntersectionGeometryStorage< LocalGeometryImpl > LocalGeometryStorageType;

    //! constructor creating an empty ALU2dGridIntersectionIterator
    ALU2dGridIntersectionBase ( const Factory &factory, int wLevel );

    //! copy constructor
    ALU2dGridIntersectionBase ( const This &other );

    //! assignment operator
    const This &operator= ( const This &other );

    //! check whether entities are the same or whether iterator is done
    bool equals ( const This &other ) const;

    //! return level of inside(entity)
    int level () const;

    //! return true if intersection is with boundary
    bool boundary() const;

    //! return boundary type
    int boundaryId () const;

    //! return the boundary segment index
    size_t boundarySegmentIndex() const;

    bool conforming () const { return current.conforming_; }

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

    int twistInSelf () const DUNE_DEPRECATED { return twistInInside(); }
    int twistInNeighbor () const DUNE_DEPRECATED { return twistInOutside(); }

    NormalType outerNormal ( const LocalCoordinate &local ) const;
    NormalType integrationOuterNormal ( const LocalCoordinate &local ) const;
    NormalType unitOuterNormal ( const LocalCoordinate &local ) const;

    NormalType centerUnitOuterNormal () const;

    LocalGeometry geometryInInside () const;
    LocalGeometry geometryInOutside () const;
    Geometry geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const;

  protected:
    const Grid &grid () const { return factory_.grid(); }

    void invalidate ();
    void done ( HElementType *inside );

    int walkLevel () const { return walkLevel_; }

    //! return true if intersection is with boundary
    void checkValid () ;

  public:
    IntersectionInfo current;

  protected:
    // the local geometries
    mutable GeometryImpl intersectionGlobal_;
    mutable LocalGeometryImpl intersectionSelfLocal_;
    mutable LocalGeometryImpl intersectionNeighborLocal_;

    // reference to factory
    const Factory &factory_;
    const LocalGeometryStorageType &localGeomStorage_;
    int walkLevel_;
  };



  // ALU2dGridIntersectionInfo
  // -------------------------

  template< class Grid >
  struct ALU2dGridIntersectionInfo
  {
    static const ALU2DSPACE ElementType elementType = Grid::elementType;

    static const int dimension = Grid::dimension;
    static const int dimensionworld  = Grid::dimensionworld;

    typedef typename ALU2dImplTraits< dimensionworld, elementType >::HElementType HElement;
    typedef typename ALU2dImplTraits< dimensionworld, elementType >::ThinelementType ThinElement;
    typedef typename ALU2dImplTraits< dimensionworld, elementType >::HBndElType HBndElement;

    explicit ALU2dGridIntersectionInfo ( HElement *inside = nullptr )
      : index_( 0 ), useOutside_( false )
    {
      setInside( inside );
      setOutside( nullptr, -1 );
    }

    HElement *inside () const { return inside_; }

    int nFaces () const { return nFaces_; }

    bool isBoundary () const
    {
      assert( inside() && (index() < nFaces()) && inside()->neighbour( index() ) );
      return inside()->neighbour( index() )->thinis( ThinElement::bndel_like );
    }

    HBndElement *boundary () const
    {
      assert( isBoundary() );
      return (HBndElement *)inside()->neighbour( index() );
    }

    HElement *outside () const { return outside_; }

    int index () const { return index_; }

    int opposite () const { return opposite_; }

    bool conforming () const { return conforming_; }

    void setInside ( HElement *inside )
    {
      inside_ = inside;
      nFaces_ = (inside ? inside->numfaces() : 0);
    }

    void setOutside ( HElement *outside, int opposite )
    {
      outside_ = outside;
      opposite_ = opposite;
    }

  private:
    HElement *inside_;
    HElement *outside_;
    int nFaces_;
    int opposite_;

  public:
    mutable int index_;
    mutable bool useOutside_;
    bool conforming_;
  };



  // ALU2dGridLeafIntersection
  // -------------------------

  template< class Grid >
  class ALU2dGridLeafIntersection
    : public ALU2dGridIntersectionBase< Grid, ALU2dGridIntersectionInfo< Grid > >
  {
    typedef ALU2dGridLeafIntersection< Grid > This;
    typedef ALU2dGridIntersectionBase< Grid, ALU2dGridIntersectionInfo< Grid > > Base;

  public:
    typedef typename Base::Factory Factory;

    ALU2dGridLeafIntersection ( const Factory &factory, int wLevel )
      : Base( factory, wLevel )
    {}
  };



  // ALU2dGridLevelIntersection
  // --------------------------

  template< class Grid >
  class ALU2dGridLevelIntersection
    : public ALU2dGridIntersectionBase< Grid, ALU2dGridIntersectionInfo< Grid > >
  {
    typedef ALU2dGridLevelIntersection< Grid > This;
    typedef ALU2dGridIntersectionBase< Grid, ALU2dGridIntersectionInfo< Grid > > Base;

  public:
    typedef typename Base::Factory Factory;

    ALU2dGridLevelIntersection ( const Factory &factory, int wLevel )
      : Base( factory, wLevel )
    {}
  };



  // ALU2dGridLevelIntersectionIterator
  // ----------------------------------

  template< class Grid >
  class ALU2dGridLevelIntersectionIterator
  {
    typedef ALU2dGridLevelIntersectionIterator< Grid > This;

    static const ALU2DSPACE ElementType eltype = Grid::elementType;

  public:
    static const int dimension = Grid::dimension;
    static const int dimensionworld  = Grid::dimensionworld;

  private:
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::ThinelementType ThinelementType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::HBndElType HBndElType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::PeriodicBndElType PeriodicBndElType;

    typedef ALU2dGridLevelIntersection< Grid > IntersectionImpl;
    typedef typename IntersectionImpl::IntersectionInfo IntersectionInfo;

    typedef typename IntersectionImpl::HElementType HElementType;
    typedef std::pair< HElementType *, int > OutsideInfo;

  public:
    //! type of the intersection
    typedef Dune::Intersection< Grid, ALU2dGridLevelIntersection > Intersection;

    typedef typename Grid::GridObjectFactoryType Factory;
    typedef ALUMemoryProvider< This > StorageType;

    typedef typename Grid::template Codim<0>::Entity Entity;

    typedef typename Grid::template Codim<1>::Geometry Geometry;
    typedef typename Grid::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dimension,Grid> EntityImp;
    typedef ALU2dGridGeometry<dimension-1,dimensionworld,Grid> GeometryImp;
    typedef ALU2dGridGeometry<dimension-1,dimension,Grid> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimensionworld> NormalType;
    typedef ALU2dGridEntityPointer<0,Grid> EntityPointer;

    /** \brief constructor
     *
     *  \param[in]  factory  factory for grid objects
     *  \param[in]  el       pointer to the inside HElement
     *  \param[in]  level    level on which we want neighbors
     *  \param[in]  end      shall this bee the end iterator?
     */
    ALU2dGridLevelIntersectionIterator ( const Factory &factory, HElementType *el, int wLevel, bool end = true );

    //! copy constructor
    ALU2dGridLevelIntersectionIterator ( const This &other );

    //! assignment operator
    const This &operator= ( const This &other );

    //! increment iterator
    void increment ();

    bool equals ( const This &other ) const { return intersectionImpl().equals( other.intersectionImpl() ); }

    const Intersection &dereference () const { return intersection_; }

  protected:
    bool isConform () const
    {
      return (!current().outside() || (current().outside() == current().inside()->neighbour( current().index() )));
    }

  private:
    void doIncrement ();

    void addNeighboursToStack ();

    static int getOppositeInFather ( int nrInChild, int nrOfChild );
    static int getOppositeInChild ( int nrInFather, int nrOfChild );

    void setupIntersection ();

  protected:
    const IntersectionImpl &intersectionImpl () const { return Grid::getRealImplementation( intersection_ ); }
    IntersectionImpl &intersectionImpl () { return Grid::getRealImplementation( intersection_ ); }

    const IntersectionInfo &current () const { return intersectionImpl().current; }
    IntersectionInfo &current () { return intersectionImpl().current; }

    int walkLevel () const { return intersectionImpl().walkLevel(); }

  private:
    Intersection intersection_;
    std::stack< OutsideInfo > nbStack_;
  };



  // ALU2dGridLeafIntersectionIterator
  // ---------------------------------

  template< class Grid >
  class ALU2dGridLeafIntersectionIterator
  {
    typedef ALU2dGridLeafIntersectionIterator< Grid > This;

    static const ALU2DSPACE ElementType eltype = Grid::elementType;

    typedef ALU2dGridLeafIntersection< Grid > IntersectionImpl;
    typedef typename IntersectionImpl::IntersectionInfo IntersectionInfo;

    typedef typename IntersectionImpl::HElementType HElementType;
    typedef std::pair< HElementType *, int > OutsideInfo;

  public:
    static const int dimension = Grid::dimension;
    static const int dimensionworld  = Grid::dimensionworld;

    typedef typename Grid::GridObjectFactoryType Factory;
    typedef ALUMemoryProvider< This > StorageType;

  private:
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::ThinelementType ThinelementType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::HBndElType HBndElType;
    typedef typename ALU2dImplTraits< dimensionworld, eltype >::PeriodicBndElType PeriodicBndElType;

  public:
    //! type of the intersection
    typedef Dune::Intersection< Grid, Dune::ALU2dGridLeafIntersection > Intersection;

    typedef typename Grid::template Codim<0>::Entity Entity;
    typedef typename Grid::template Codim<1>::Geometry Geometry;
    typedef typename Grid::template Codim<1>::LocalGeometry LocalGeometry;
    typedef ALU2dGridEntity<0,dimension,Grid> EntityImp;
    typedef ALU2dGridGeometry<dimension-1,dimensionworld,Grid> GeometryImp;
    typedef ALU2dGridGeometry<dimension-1,dimension,Grid> LocalGeometryImp;
    typedef FieldVector<alu2d_ctype, dimensionworld> NormalType;
    typedef ALU2dGridEntityPointer<0,Grid> EntityPointer;

    /** \brief constructor
     *
     *  \param[in]  factory  factory for grid objects
     *  \param[in]  el       pointer to the inside HElement
     *  \param[in]  level    level on which we want neighbors
     *  \param[in]  end      shall this bee the end iterator?
     */
    ALU2dGridLeafIntersectionIterator ( const Factory &factory, HElementType *el, int wLevel, bool end = true );

    //! copy constructor
    ALU2dGridLeafIntersectionIterator ( const This &other );

    //! assignment operator
    const This &operator= ( const This &other );

    //! increment iterator
    void increment ();

    bool equals ( const This &other ) const { return intersectionImpl().equals( other.intersectionImpl() ); }

    const Intersection &dereference () const { return intersection_; }

  private:
    void doIncrement ();

    void setupIntersection ();

  protected:
    const IntersectionImpl &intersectionImpl () const { return Grid::getRealImplementation( intersection_ ); }
    IntersectionImpl &intersectionImpl () { return Grid::getRealImplementation( intersection_ ); }

    const IntersectionInfo &current () const { return intersectionImpl().current; }
    IntersectionInfo &current () { return intersectionImpl().current; }

    int walkLevel () const { return intersectionImpl().walkLevel(); }

  private:
    Intersection intersection_;
    std::stack< OutsideInfo > nbStack_;
  };

} // namespace Dune

#include "intersection_imp.cc"
#if COMPILE_ALU2DGRID_INLINE
  #include "intersection.cc"
#endif

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_HH
