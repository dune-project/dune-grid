// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTION_HH
#define DUNE_GEOGRID_INTERSECTION_HH

#include <dune/grid/geometrygrid/entitypointer.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  class GeometryGrid;



  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid >
    class LeafIntersection;

    template< class Grid >
    class LevelIntersection;



    // Intersection
    // ------------

    template< class Grid, class HostIntersection >
    class Intersection
    {
      typedef typename HostIntersection::Geometry HostGeometry;
      typedef typename HostIntersection::LocalGeometry HostLocalGeometry;

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;

    private:
      typedef typename GenericGeometry::GlobalGeometryTraits< Grid >::IntersectionCoordVector CoordVector;

      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry::ImplementationType GeometryImpl;

      typedef typename MakeableInterfaceObject< ElementGeometry >::ImplementationType ElementGeometryImpl;

    public:
      Intersection ( const Grid &grid, const ElementGeometry &insideGeo )
        : grid_( &grid ),
          insideGeo_( Grid::getRealImplementation( insideGeo ) ),
          hostIntersection_( 0 ),
          geo_( GeometryImpl( grid.allocator() ) )
      {}

      Intersection ( const Intersection &other )
        : grid_( other.grid_ ),
          insideGeo_( Grid::getRealImplementation( other.insideGeo_ ) ),
          hostIntersection_( 0 ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      const Intersection &operator= ( const Intersection &other )
      {
        grid_ = other.grid_;
        Grid::getRealImplementation( insideGeo_ ) = Grid::getRealImplementation( other.insideGeo_ );
        invalidate();
        return *this;
      }

      operator bool () const { return bool( hostIntersection_ ); }

      EntityPointer inside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().inside() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().outside() );
      }

      bool boundary () const { return hostIntersection().boundary(); }

      bool conforming () const { return hostIntersection().conforming(); }

      bool neighbor () const { return hostIntersection().neighbor(); }

      int boundaryId () const { return hostIntersection().boundaryId(); }

      size_t boundarySegmentIndex () const
      {
        return hostIntersection().boundarySegmentIndex();
      }

      const LocalGeometry &geometryInInside () const
      {
        return hostIntersection().geometryInInside();
      }

      const LocalGeometry &geometryInOutside () const
      {
        return hostIntersection().geometryInOutside();
      }

      const Geometry &geometry () const
      {
        GeometryImpl &geo = Grid::getRealImplementation( geo_ );
        if( !geo )
        {
          const LocalGeometry &localGeo = geometryInInside();
          CoordVector coords( insideGeometry(), localGeo );
          geo = GeometryImpl( type(), coords, grid().allocator() );
        }
        return geo_;
      }

      GeometryType type () const { return hostIntersection().type(); }

      unsigned int topologyId () const DUNE_DEPRECATED
      {
        return type().id();
      }

      int indexInInside () const
      {
        return hostIntersection().indexInInside();
      }

      int indexInOutside () const
      {
        return hostIntersection().indexInOutside();
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        const ElementGeometryImpl &geo = Grid::getRealImplementation( insideGeometry() );

        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( geo.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        const typename ElementGeometryImpl::JacobianInverseTransposed &jit = geo.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        normal *= ctype( 1 ) / jit.det();
        //normal *= geo.integrationElement( x );
        return normal;
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return integrationOuterNormal( local );
      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        FieldVector< ctype, dimensionworld > normal = outerNormal( local );
        normal *= (ctype( 1 ) / normal.two_norm());
        return normal;
      }

      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        const GenericReferenceElement< ctype, dimension-1 > &refFace
          = GenericReferenceElements< ctype, dimension-1 >::general( type() );
        return unitOuterNormal( refFace.position( 0, 0 ) );
      }

      const HostIntersection &hostIntersection () const
      {
        assert( *this );
        return *hostIntersection_;
      }

      const Grid &grid () const { return *grid_; }

      void invalidate ()
      {
        hostIntersection_ = 0;
        Grid::getRealImplementation( geo_ ) = GeometryImpl( grid().allocator() );
      }

      void initialize ( const HostIntersection &hostIntersection )
      {
        assert( !(*this) );
        hostIntersection_ = &hostIntersection;
      }

    private:
      const ElementGeometry &insideGeometry () const { return insideGeo_; }

      const Grid *grid_;
      ElementGeometry insideGeo_;
      const HostIntersection *hostIntersection_;
      mutable MakeableGeometry geo_;
    };



    // LeafIntersection
    // ----------------

    template< class HostGrid, class CoordFunction, class Allocator >
    class LeafIntersection< const GeometryGrid< HostGrid, CoordFunction, Allocator > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction, Allocator >,
            typename HostGrid::Traits::LeafIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;
      typedef typename HostGrid::Traits::LeafIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

    public:
      typedef typename Base::ElementGeometry ElementGeometry;

      LeafIntersection ( const Grid &grid, const ElementGeometry &insideGeo )
        : Base( grid, insideGeo )
      {}
    };



    // LevelIntersection
    // -----------------

    template< class HostGrid, class CoordFunction, class Allocator >
    class LevelIntersection< const GeometryGrid< HostGrid, CoordFunction, Allocator > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction, Allocator >,
            typename HostGrid::Traits::LevelIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;
      typedef typename HostGrid::Traits::LevelIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

    public:
      typedef typename Base::ElementGeometry ElementGeometry;

      LevelIntersection ( const Grid &grid, const ElementGeometry &insideGeo )
        : Base( grid, insideGeo )
      {}
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_INTERSECTION_HH
