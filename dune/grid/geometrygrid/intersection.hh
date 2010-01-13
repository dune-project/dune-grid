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

  template< class HostGrid, class CoordFunction >
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

    template< class HostGrid, class CoordFunction, class HostIntersection >
    class Intersection< const GeometryGrid< HostGrid, CoordFunction >, HostIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostIntersection :: Geometry HostGeometry;
      typedef typename HostIntersection :: LocalGeometry HostLocalGeometry;

    public:
      typedef typename Grid :: ctype ctype;

      enum { dimension = Grid :: dimension };
      enum { dimensionworld = Grid :: dimensionworld };

      typedef typename Grid :: template Codim< 0 > :: Entity Entity;
      typedef typename Grid :: template Codim< 0 > :: EntityPointer EntityPointer;
      typedef typename Grid :: template Codim< 1 > :: Geometry Geometry;
      typedef typename Grid :: template Codim< 1 > :: LocalGeometry LocalGeometry;

    private:
      typedef typename GenericGeometry :: GlobalGeometryTraits<Grid> :: IntersectionCoordVector CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry :: ImplementationType GeometryImpl;

      const EntityPointer *inside_;
      const HostIntersection *hostIntersection_;
      mutable MakeableGeometry geo_;

    public:
      Intersection ()
        : geo_( GeometryImpl() )
      {}

      Intersection ( const Intersection &other )
        : inside_( other.inside_ ),
          hostIntersection_( other.hostIntersection_ ),
          geo_( GeometryImpl() )
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
          geo = GeometryImpl( localGeo.type(), coords );
        }
        return geo_;
      }

      GeometryType type () const
      {
        return hostIntersection().type();
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
        typedef typename Grid :: template Codim< 0 > :: Geometry Geometry;
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

      const HostIntersection &hostIntersection () const
      {
        assert( isValid() );
        return *hostIntersection_;
      }

      void initialize( const EntityPointer &inside, const HostIntersection &hostIntersection )
      {
        inside_ = &inside;
        hostIntersection_ = &hostIntersection;
        Grid :: getRealImplementation( geo_ ) = GeometryImpl();
      }

    protected:
      const Grid &grid () const
      {
        return Grid :: getRealImplementation( inside() ).grid();
      }

      bool isValid () const
      {
        return (hostIntersection_ != 0);
      }

      void invalidate ()
      {
        hostIntersection_ = 0;
      }
    };


    // LeafIntersection
    // ----------------

    template< class HostGrid, class CoordFunction >
    class LeafIntersection< const GeometryGrid< HostGrid, CoordFunction > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction >,
            typename HostGrid :: Traits :: LeafIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;
      typedef typename HostGrid :: Traits :: LeafIntersection HostIntersection;

      template< class > friend class IntersectionWrapper;
    };



    // LevelIntersection
    // -----------------

    template< class HostGrid, class CoordFunction >
    class LevelIntersection< const GeometryGrid< HostGrid, CoordFunction > >
      : public Intersection
        < const GeometryGrid< HostGrid, CoordFunction >,
            typename HostGrid :: Traits :: LevelIntersection >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;
      typedef typename HostGrid :: Traits :: LevelIntersection HostIntersection;

      template< class > friend class IntersectionWrapper;
    };



    // IntersectionWrapper
    // -------------------

    template< class Intersection >
    class IntersectionWrapper
      : public Intersection
    {
      typedef Intersection Base;

    protected:
      using Base :: getRealImp;

    public:
      typedef typename Intersection :: ImplementationType Implementation;

      typedef typename Implementation :: EntityPointer EntityPointer;
      typedef typename Implementation :: HostIntersection HostIntersection;

      IntersectionWrapper ()
        : Base( Implementation() )
      {}

      void initialize( const EntityPointer &inside,
                       const HostIntersection &hostIntersection )
      {
        getRealImp().initialize( inside, hostIntersection );
      }
    };

  }

}

#endif
