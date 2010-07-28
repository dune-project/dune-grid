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

  template< class HostGrid, class CoordFunction, class Numbering, class Allocator >
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

      typedef typename Traits::IntersectionNumbering Numbering;

    public:
      Intersection ( const Grid &grid, const ElementGeometry &insideGeo, const HostIntersection &hostIntersection )
        : grid_( &grid ),
          insideGeo_( Grid::getRealImplementation( insideGeo ) ),
          hostIntersection_( &hostIntersection ),
          numbering_( grid.numbering()[ hostIntersection ] ),
          geo_( GeometryImpl( grid.allocator() ) )
      {}

      Intersection ( const Intersection &other )
        : grid_( other.grid_ ),
          insideGeo_( Grid::getRealImplementation( other.insideGeo_ ) ),
          hostIntersection_( other.hostIntersection_ ),
          numbering_( other.numbering_ ),
          geo_( Grid::getRealImplementation( other.geo_ ) )
      {}

      EntityPointer inside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().inside() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().outside() );
      }

      bool boundary () const
      {
        return hostIntersection().boundary();
      }

      bool conforming () const
      {
        return hostIntersection().conforming();
      }

      bool neighbor () const
      {
        return hostIntersection().neighbor();
      }

      int boundaryId () const
      {
        return hostIntersection().boundaryId();
      }

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
          geo = GeometryImpl( topologyId(), coords, grid().allocator() );
        }
        return geo_;
      }

      GeometryType type () const
      {
        return hostIntersection().type();
      }

      unsigned int topologyId () const
      {
        return GenericGeometry::topologyId( type() );
      }

      int indexInInside () const
      {
        const int i = hostIntersection().indexInInside();
        return numbering_.template map< Numbering::Inside >( i );
      }

      int indexInOutside () const
      {
        const int i = hostIntersection().indexInOutside();
        return numbering_.template map< Numbering::Outside >( i );
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        typedef typename Grid::template Codim< 0 >::Geometry Geometry;
        const Geometry &geo = insideGeometry();

        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension>::general( geo.type() );

        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        const typename Geometry::Jacobian &jit = geo.jacobianInverseTransposed( x );
        const FieldVector< ctype, dimension > &refNormal = refElement.volumeOuterNormal( indexInInside() );

        FieldVector< ctype, dimensionworld > normal;
        jit.mv( refNormal, normal );
        normal *= geo.integrationElement( x );
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
        assert( isValid() );
        return *hostIntersection_;
      }

    protected:
      const Grid &grid () const
      {
        return *grid_;
      }

      bool isValid () const
      {
        return (hostIntersection_ != 0);
      }

      void invalidate ()
      {
        hostIntersection_ = 0;
      }

    private:
      const ElementGeometry &insideGeometry () const
      {
        return insideGeo_;
      }

      const Grid *grid_;
      ElementGeometry insideGeo_;
      const HostIntersection *hostIntersection_;
      Numbering numbering_;
      mutable MakeableGeometry geo_;
    };



    // LeafIntersection
    // ----------------

    template< class HostGrid, class CoordFunction, class Numbering, class Allocator >
    class LeafIntersection< const GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator >,
            typename HostGrid::Traits::LeafIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > Grid;
      typedef typename HostGrid::Traits::LeafIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

    public:
      typedef typename Base::ElementGeometry ElementGeometry;

      LeafIntersection ( const Grid &grid, const ElementGeometry &insideGeo, const HostIntersection &hostIntersection )
        : Base( grid, insideGeo, hostIntersection )
      {}
    };



    // LevelIntersection
    // -----------------

    template< class HostGrid, class CoordFunction, class Numbering, class Allocator >
    class LevelIntersection< const GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator >,
            typename HostGrid::Traits::LevelIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > Grid;
      typedef typename HostGrid::Traits::LevelIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

    public:
      typedef typename Base::ElementGeometry ElementGeometry;

      LevelIntersection ( const Grid &grid, const ElementGeometry &insideGeo, const HostIntersection &hostIntersection )
        : Base( grid, insideGeo, hostIntersection )
      {}
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_INTERSECTION_HH
