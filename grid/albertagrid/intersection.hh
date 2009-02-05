// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_HH
#define DUNE_ALBERTA_INTERSECTION_HH

#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dune/grid/albertagrid/transformation.hh>
#include <dune/grid/albertagrid/agmemory.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/geometry.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, int dim, class GridImp >
  class AlbertaGridEntity;



  // AlbertaGridIntersectionIterator
  // -------------------------------

  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1
     These neighbors are accessed via a IntersectionIterator. This allows the implementation of
     non-matching meshes. The number of neigbors may be different from the number of faces
     of an element!
   */
  template< class GridImp >
  class AlbertaGridIntersectionIterator
  {
    typedef AlbertaGridIntersectionIterator This;

    friend class AlbertaGridEntity< 0, GridImp::dimension, GridImp >;

  public:
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! know your own dimension
    static const int dimension = GridImp::dimension;
    //! know your own dimension of world
    static const int dimensionworld = GridImp::dimensionworld;

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    typedef FieldVector< ctype, GridImp::dimensionworld > NormalVector;
    typedef FieldVector< ctype, GridImp::dimension-1 > LocalCoordType;

    typedef Dune::Intersection< GridImp, Dune::AlbertaGridIntersectionIterator >
    Intersection;
    typedef This ImplementationType;

    typedef AGMemoryProvider< This > StorageType;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  private:
    typedef AlbertaGridEntity< 0, dimension, GridImp > EntityImp;
    typedef AlbertaGridGeometry< dimension-1, dimensionworld, GridImp > GeometryImp;
    typedef AlbertaGridGeometry< dimension-1, dimension, GridImp > LocalGeometryImp;

    struct LocalCoordReader;

  public:
    //! The default Constructor
    AlbertaGridIntersectionIterator ( const GridImp &grid, int level );

    //! The copy constructor
    AlbertaGridIntersectionIterator( const This &other );

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }

    //! equality
    bool equals ( const This &other ) const;

    //! increment
    void increment();

    //! access neighbor
    EntityPointer outside() const;

    //! access element where IntersectionIterator started
    EntityPointer inside() const;

    //! assignment operator, implemented because default does not the right thing
    void assign ( const This &other );

    //! return true if intersection is with boundary.
    bool boundary () const;

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const;

    //! return information about the Boundary
    int boundaryId () const;

    //! return true if intersection is conform.
    bool conforming () const;

    AlbertaTransformation transformation () const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &intersectionSelfLocal () const;

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of neighbor
     */
    const LocalGeometry &intersectionNeighborLocal () const;

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in GLOBAL coordinates of the element where iteration started.
     */
    const Geometry &intersectionGlobal () const;

    GeometryType type () const;

    //! local number of codim 1 entity in self where intersection is contained in
    int numberInSelf () const;

    //! local number of codim 1 entity in neighbor where intersection is contained in
    int numberInNeighbor () const;

    //! twist of the face seen from the inner element
    int twistInSelf () const;

    //! twist of the face seen from the outer element
    int twistInNeighbor () const;

    const NormalVector unitOuterNormal ( const LocalCoordType &local ) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalVector outerNormal ( const LocalCoordType &local ) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalVector integrationOuterNormal ( const LocalCoordType &local ) const;

    //**********************************************************
    //  private methods
    //**********************************************************

    // reset IntersectionIterator
    void first ( const EntityImp &entity, int level );

    // calls EntityPointer done and sets done_ to true
    void done ();

  private:
    // returns true if actual neighbor has same level
    bool neighborHasSameLevel () const;

    //! setup the virtual neighbor
    void setupVirtEn () const;

    //! calculate normal to current face
    void calcOuterNormal ( NormalVector &n ) const;

    ////////////////////////////////////////////////
    // private member variables
    ////////////////////////////////////////////////

    //! know the grid were im coming from
    const GridImp &grid_;

    //! count on which neighbor we are lookin' at
    mutable int neighborCount_;

    ElementInfo elementInfo_;

    // the objects holding the real implementations
    mutable MakeableInterfaceObject< LocalGeometry > fakeNeighObj_;
    mutable MakeableInterfaceObject< LocalGeometry > fakeSelfObj_;
    mutable MakeableInterfaceObject< Geometry > neighGlobObj_;

    //! ElementInfo to store the information of the neighbor if needed
    mutable ElementInfo neighborInfo_;
  };

}

#endif
