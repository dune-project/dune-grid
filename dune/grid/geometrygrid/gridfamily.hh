// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRIDFAMILY_HH
#define DUNE_GEOGRID_GRIDFAMILY_HH

#include <dune/common/static_assert.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/entity.hh>
#include <dune/grid/geometrygrid/entityseed.hh>
#include <dune/grid/geometrygrid/entitypointer.hh>
#include <dune/grid/geometrygrid/geometry.hh>
#include <dune/grid/geometrygrid/gridview.hh>
#include <dune/grid/geometrygrid/intersection.hh>
#include <dune/grid/geometrygrid/intersectioniterator.hh>
#include <dune/grid/geometrygrid/iterator.hh>
#include <dune/grid/geometrygrid/idset.hh>
#include <dune/grid/geometrygrid/indexsets.hh>

namespace Dune
{

  /** \brief namespace containing the implementations of GeometryGrid
   *  \ingroup GeoGrid
   */
  namespace GeoGrid
  {

    // ExportParams
    // ------------

    template< class HG, class CF >
    class ExportParams
    {
      static const bool isCoordFunction = isCoordFunctionInterface< typename CF::Interface >::value;
      dune_static_assert( isCoordFunction, "Invalid CoordFunction." );

    public:
      typedef HG HostGrid;
      typedef CF CoordFunction;
    };



    // GridFamily
    // ----------

    template< class HG, class CF, class Allocator >
    struct GridFamily
    {
      struct Traits
      {
        typedef GeometryGrid< HG, CF, Allocator > Grid;

        typedef HG HostGrid;
        typedef CF CoordFunction;

        typedef typename HostGrid::ctype ctype;

        static const int dimension = HostGrid::dimension;
        static const int dimensionworld = CoordFunction::dimRange;

        typedef Dune::Intersection< const Grid, GeoGrid::Intersection< const Grid, typename HostGrid::LeafIntersection > > LeafIntersection;
        typedef Dune::Intersection< const Grid, GeoGrid::Intersection< const Grid, typename HostGrid::LevelIntersection > > LevelIntersection;

        typedef Dune::IntersectionIterator
        < const Grid, GeoGrid::IntersectionIterator< const Grid, typename HostGrid::LeafIntersectionIterator >, GeoGrid::Intersection< const Grid, typename HostGrid::LeafIntersection > >
        LeafIntersectionIterator;
        typedef Dune::IntersectionIterator
        < const Grid, GeoGrid::IntersectionIterator< const Grid, typename HostGrid::LevelIntersectionIterator >, GeoGrid::Intersection< const Grid, typename HostGrid::LevelIntersection > >
        LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const Grid, GeoGrid::HierarchicIterator< const Grid > >
        HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune::GeoGrid::Geometry< dimension-codim, dimensionworld, const Grid > GeometryImpl;
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, Dune::GeoGrid::Geometry > Geometry;
          typedef typename HostGrid::template Codim< codim >::LocalGeometry LocalGeometry;

          typedef GeoGrid::EntityPointerTraits< codim, const Grid > EntityPointerTraits;
          typedef GeoGrid::EntityPointer< EntityPointerTraits > EntityPointerImpl;
          typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;
          typedef typename EntityPointerTraits::Entity Entity;

          typedef Dune::EntitySeed< const Grid, GeoGrid::EntitySeed< codim, const Grid > > EntitySeed;

          template< PartitionIteratorType pitype >
          struct Partition
          {
            typedef GeoGrid::IteratorTraits< typename HostGrid::LeafGridView, codim, pitype, const Grid > LeafIteratorTraits;
            typedef Dune::EntityIterator< codim, const Grid, GeoGrid::Iterator< LeafIteratorTraits > > LeafIterator;

            typedef GeoGrid::IteratorTraits< typename HostGrid::LevelGridView, codim, pitype, const Grid > LevelIteratorTraits;
            typedef Dune::EntityIterator< codim, const Grid, GeoGrid::Iterator< LevelIteratorTraits > > LevelIterator;
          };

          typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
          typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        };

        typedef GeoGrid::IndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet > LeafIndexSet;
        typedef GeoGrid::IndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet > LevelIndexSet;

        typedef GeoGrid::IdSet< const Grid, typename HostGrid::Traits::GlobalIdSet >
        GlobalIdSet;
        typedef GeoGrid::IdSet< const Grid, typename HostGrid::Traits::LocalIdSet >
        LocalIdSet;

        typedef typename HostGrid::Traits::CollectiveCommunication CollectiveCommunication;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::GridView< GeoGrid::GridViewTraits< typename HostGrid::LeafGridView, CoordFunction, Allocator, pitype > >
          LeafGridView;
          typedef Dune::GridView< GeoGrid::GridViewTraits< typename HostGrid::LevelGridView, CoordFunction, Allocator, pitype > >
          LevelGridView;
        };
      };
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GRIDFAMILY_HH
