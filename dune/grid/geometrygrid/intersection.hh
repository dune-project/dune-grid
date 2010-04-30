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

    template< class Grid, class HostIntersection >
    class Intersection;

    template< class Grid >
    class LeafIntersection;

    template< class Grid >
    class LevelIntersection;

    template< class Intersection >
    class IntersectionWrapper;



    // Intersection
    // ------------

    template< class HostGrid, class CoordFunction, class Numbering, class Allocator, class HostIntersection >
    class Intersection< const GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator >, HostIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > Grid;

      typedef typename HostIntersection::Geometry HostGeometry;
      typedef typename HostIntersection::LocalGeometry HostLocalGeometry;

    public:
      typedef typename Grid::ctype ctype;

      enum { dimension = Grid::dimension };
      enum { dimensionworld = Grid::dimensionworld };

      typedef typename Grid::template Codim< 0 >::Entity Entity;
      typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Grid::template Codim< 1 >::Geometry Geometry;
      typedef typename Grid::template Codim< 1 >::LocalGeometry LocalGeometry;

    private:
      typedef typename GenericGeometry::GlobalGeometryTraits< Grid >::IntersectionCoordVector CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry::ImplementationType GeometryImpl;

      typedef typename Grid::Traits::IntersectionNumbering IntersectionNumbering;

    public:
      Intersection ( const EntityPointer &inside, const HostIntersection &hostIntersection )
        : inside_( &inside ),
          hostIntersection_( &hostIntersection ),
          numbering_( grid().numbering()[ hostIntersection ] ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      Intersection ( const Intersection &other )
        : inside_( other.inside_ ),
          hostIntersection_( other.hostIntersection_ ),
          numbering_( other.numbering_ ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      const EntityPointer &inside () const
      {
        return *inside_;
      }

      EntityPointer outside () const
      {
        typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
        typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
        return MakeableEntityPointer( EntityPointerImpl( grid(), hostIntersection().outside() ) );
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
          CoordVector coords( inside()->geometry(), localGeo );
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
        return hostIntersection().indexInInside();
      }

      int indexInOutside () const
      {
        return hostIntersection().indexInOutside();
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        typedef typename Grid::template Codim< 0 >::Geometry Geometry;
        const Geometry &geo = inside()->geometry();
        FieldVector< ctype, dimension > x( geometryInInside().global( local ) );
        return Grid::getRealImplementation( geo ).normal( indexInInside(), x );
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
        return Grid::getRealImplementation( inside() ).grid();
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
      const EntityPointer *inside_;
      const HostIntersection *hostIntersection_;
      IntersectionNumbering numbering_;
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
      template< class > friend class IntersectionWrapper;

      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > Grid;
      typedef typename HostGrid::Traits::LeafIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

      typedef typename Base::EntityPointer EntityPointer;

    public:
      LeafIntersection ( const EntityPointer &inside,
                         const HostIntersection &hostIntersection )
        : Base( inside, hostIntersection )
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
      template< class > friend class IntersectionWrapper;

      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Allocator > Grid;
      typedef typename HostGrid::Traits::LevelIntersection HostIntersection;

      typedef Intersection< const Grid, HostIntersection > Base;

      typedef typename Base::EntityPointer EntityPointer;

    public:
      LevelIntersection ( const EntityPointer &inside,
                          const HostIntersection &hostIntersection )
        : Base( inside, hostIntersection )
      {}
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_INTERSECTION_HH
