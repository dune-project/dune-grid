// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_HH
#define DUNE_ALBERTA_INTERSECTION_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/grid/albertagrid/transformation.hh>
#include <dune/grid/albertagrid/agmemory.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/geometry.hh>

#define ALBERTA_CACHED_LOCAL_INTERSECTION_GEOMETRIES 1

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, int dim, class GridImp >
  class AlbertaGridEntity;



  // AlbertaGridIntersectionBase
  // ---------------------------

  template< class Grid >
  class AlbertaGridIntersectionBase
    : public SmallObject
  {
    typedef AlbertaGridIntersectionBase< Grid > This;

  public:
    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

    typedef FieldVector< ctype, dimensionworld > NormalVector;
    typedef FieldVector< ctype, dimension-1 > LocalCoordType;

    typedef typename Grid::template Codim< 0 >::Entity Entity;
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    typedef typename Grid::template Codim< 1 >::Geometry Geometry;
    typedef typename Grid::template Codim< 1 >::LocalGeometry LocalGeometry;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  protected:
    typedef AlbertaGridEntity< 0, dimension, Grid > EntityImp;
    typedef AlbertaGridGeometry< dimension-1, dimensionworld, Grid > GeometryImp;
    typedef AlbertaGridGeometry< dimension-1, dimension, Grid > LocalGeometryImp;

  public:
    AlbertaGridIntersectionBase ( const EntityImp &entity, const int oppVertex );

    EntityPointer inside () const;

    bool boundary () const;
    int boundaryId () const;

    int indexInInside () const;

    GeometryType type () const;

    const NormalVector integrationOuterNormal ( const LocalCoordType &local ) const;
    const NormalVector outerNormal ( const LocalCoordType &local ) const;
    const NormalVector unitOuterNormal ( const LocalCoordType &local ) const;


    const Grid &grid () const;
    const ElementInfo &elementInfo () const;

  protected:
    const Grid *grid_;
    ElementInfo elementInfo_;
    int oppVertex_;
  };



  // AlbertaGridLeafIntersection
  // ---------------------------

  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1
     These neighbors are accessed via a IntersectionIterator. This allows the implementation of
     non-matching meshes. The number of neigbors may be different from the number of faces
     of an element!
   */
  template< class GridImp >
  class AlbertaGridLeafIntersection
    : public AlbertaGridIntersectionBase< GridImp >
  {
    typedef AlbertaGridLeafIntersection< GridImp > This;
    typedef AlbertaGridIntersectionBase< GridImp > Base;

    friend class AlbertaGridEntity< 0, GridImp::dimension, GridImp >;

  public:
    typedef This ImplementationType;

    static const int dimension = Base::dimension;

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    typedef typename Base::NormalVector NormalVector;
    typedef typename Base::LocalCoordType LocalCoordType;

    typedef typename Base::Entity Entity;
    typedef typename Base::EntityPointer EntityPointer;

    typedef typename Base::Geometry Geometry;
    typedef typename Base::LocalGeometry LocalGeometry;

    typedef typename Base::ElementInfo ElementInfo;

  protected:
    typedef typename Base::EntityImp EntityImp;
    typedef typename Base::GeometryImp GeometryImp;
    typedef typename Base::LocalGeometryImp LocalGeometryImp;

    using Base::grid;
    using Base::elementInfo;

  private:
    struct GlobalCoordReader;
    struct LocalCoordReader;

  public:
    AlbertaGridLeafIntersection ( const EntityImp &entity, const int n );

    AlbertaGridLeafIntersection ( const This &other );

    This &operator= ( const This &other )
    {
      assign( other );
      return *this;
    }

    bool equals ( const This &other ) const;

    void next ();

    EntityPointer outside () const;

    void assign ( const This &other );

    bool neighbor () const;

    bool conforming () const;

    AlbertaTransformation transformation () const;

    const LocalGeometry &geometryInInside () const;
    const LocalGeometry &geometryInOutside () const;

    const Geometry &geometry () const;

    int indexInOutside () const;

    int twistInSelf () const;
    int twistInNeighbor () const;

  private:
    void setupVirtEn () const;

  protected:
    using Base::oppVertex_;

  private:
#if not ALBERTA_CACHED_LOCAL_INTERSECTION_GEOMETRIES
    mutable MakeableInterfaceObject< LocalGeometry > fakeNeighObj_;
    mutable MakeableInterfaceObject< LocalGeometry > fakeSelfObj_;
#endif
    mutable MakeableInterfaceObject< Geometry > neighGlobObj_;
    mutable ElementInfo neighborInfo_;
  };

}

#endif // #ifndef DUNE_ALBERTA_INTERSECTION_HH
