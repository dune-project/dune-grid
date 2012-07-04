// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INTERSECTION_HH
#define DUNE_GEOGRID_INTERSECTION_HH

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/entitypointer.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

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
      typedef GeoGrid::IntersectionCoordVector< Grid > CoordVector;

      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;

      typedef typename Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;

    public:
      explicit Intersection ( const ElementGeometry &insideGeo )
        : insideGeo_( Grid::getRealImplementation( insideGeo ) ),
          hostIntersection_( 0 ),
          geo_( grid() )
      {}

      Intersection ( const Intersection &other )
        : insideGeo_( other.insideGeo_ ),
          hostIntersection_( 0 ),
          geo_( grid() )
      {}

      const Intersection &operator= ( const Intersection &other )
      {
        insideGeo_ = other.insideGeo_;
        invalidate();
        return *this;
      }

      operator bool () const { return bool( hostIntersection_ ); }

      EntityPointer inside () const
      {
        return EntityPointerImpl( insideGeo_, hostIntersection().inside() );
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

      LocalGeometry geometryInInside () const
      {
        return hostIntersection().geometryInInside();
      }

      LocalGeometry geometryInOutside () const
      {
        return hostIntersection().geometryInOutside();
      }

      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVector coords( insideGeo_, geometryInInside() );
          geo_ = GeometryImpl( grid(), type(), coords );
        }
        return Geometry( geo_ );
      }

      GeometryType type () const { return hostIntersection().type(); }

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
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( insideGeo_.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        const typename ElementGeometryImpl::JacobianInverseTransposed &jit = insideGeo_.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        normal *= ctype( 1 ) / jit.det();
        //normal *= insideGeo_.integrationElement( x );
        return normal;
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( insideGeo_.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        const typename ElementGeometryImpl::JacobianInverseTransposed &jit = insideGeo_.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        return normal;
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

      const Grid &grid () const { return insideGeo_.grid(); }

      void invalidate ()
      {
        hostIntersection_ = 0;
        geo_ = GeometryImpl( grid() );
      }

      void initialize ( const HostIntersection &hostIntersection )
      {
        assert( !(*this) );
        hostIntersection_ = &hostIntersection;
      }

    private:
      ElementGeometryImpl insideGeo_;
      const HostIntersection *hostIntersection_;
      mutable GeometryImpl geo_;
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

      explicit LeafIntersection ( const ElementGeometry &insideGeo )
        : Base( insideGeo )
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

      explicit LevelIntersection ( const ElementGeometry &insideGeo )
        : Base( insideGeo )
      {}
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_INTERSECTION_HH
