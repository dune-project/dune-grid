// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_HH
#define DUNE_ALBERTA_INTERSECTION_HH

#include <dune/grid/common/intersection.hh>

#include <dune/grid/albertagrid/transformation.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/geometry.hh>

#if HAVE_ALBERTA

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
  {
    typedef AlbertaGridIntersectionBase< Grid > This;

  public:
    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

    typedef FieldVector< ctype, dimensionworld > NormalVector;
    typedef FieldVector< ctype, dimension-1 > LocalCoordType;

    typedef typename Grid::template Codim< 0 >::Entity Entity;

    typedef typename Grid::template Codim< 1 >::Geometry Geometry;
    typedef typename Grid::template Codim< 1 >::LocalGeometry LocalGeometry;

    typedef Alberta::ElementInfo< dimension > ElementInfo;

  protected:
    typedef AlbertaGridEntity< 0, dimension, Grid > EntityImp;

    typedef typename Grid::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename Grid::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

    struct GlobalCoordReader;
    struct LocalCoordReader;

  public:

    AlbertaGridIntersectionBase ();

    AlbertaGridIntersectionBase ( const EntityImp &entity, const int oppVertex );

    Entity inside () const;

    bool boundary () const;
    int boundaryId () const;
    size_t boundarySegmentIndex () const;

    int indexInInside () const;

    GeometryType type () const;

    NormalVector centerIntegrationOuterNormal () const;
    NormalVector centerOuterNormal () const;
    NormalVector centerUnitOuterNormal () const;

    NormalVector integrationOuterNormal ( [[maybe_unused]] const LocalCoordType &local ) const;
    NormalVector outerNormal ( [[maybe_unused]] const LocalCoordType &local ) const;
    NormalVector unitOuterNormal ( [[maybe_unused]] const LocalCoordType &local ) const;


    AlbertaTransformation transformation () const;


    const Grid &grid () const;
    const ElementInfo &elementInfo () const;

  protected:
    const Grid *grid_;
    ElementInfo elementInfo_;
    int oppVertex_;
  };



  // AlbertaGridLeafIntersection
  // ---------------------------

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
    static const int dimensionworld = Base::dimensionworld;

    typedef typename Base::NormalVector NormalVector;
    typedef typename Base::LocalCoordType LocalCoordType;

    typedef typename Base::Entity Entity;

    typedef typename Base::Geometry Geometry;
    typedef typename Base::LocalGeometry LocalGeometry;

    typedef typename Base::ElementInfo ElementInfo;

  protected:
    typedef typename Base::EntityImp EntityImp;

    typedef typename Base::GeometryImpl GeometryImpl;
    typedef typename Base::LocalGeometryImpl LocalGeometryImpl;

    typedef typename Base::GlobalCoordReader GlobalCoordReader;
    typedef typename Base::LocalCoordReader LocalCoordReader;

  public:
    using Base::grid;
    using Base::elementInfo;

    using Base::inside;

    AlbertaGridLeafIntersection () = default;

    AlbertaGridLeafIntersection ( const EntityImp &entity, int n ) : Base( entity, n ) {}

    AlbertaGridLeafIntersection ( const This &other ) : Base( other ) {}

    This &operator= ( const This &other )
    {
      *static_cast< Base * >( this ) = other;
      neighborInfo_ = ElementInfo();
      return *this;
    }

    bool operator== ( const This &other ) const { return (oppVertex_ == other.oppVertex_) && (elementInfo_ == other.elementInfo_); }
    bool operator!= ( const This &other ) const { return (oppVertex_ != other.oppVertex_) || (elementInfo_ != other.elementInfo_); }

    bool equals ( const AlbertaGridLeafIntersection& other ) const { return (*this) == other; }

    void next ();

    typename GridImp::template Codim< 0 >::Entity outside () const;

    bool neighbor () const;

    bool conforming () const { return true; }

    LocalGeometry geometryInInside () const;
    LocalGeometry geometryInOutside () const;

    Geometry geometry () const;

    int indexInOutside () const;

    int twistInInside () const { return elementInfo().template twist< 1 >( oppVertex_ ); }
    int twistInOutside () const { return elementInfo().twistInNeighbor( oppVertex_ ); }

  protected:
    using Base::elementInfo_;
    using Base::oppVertex_;

  private:
    mutable ElementInfo neighborInfo_;
  };

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_INTERSECTION_HH
