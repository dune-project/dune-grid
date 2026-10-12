// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKGEOMETRY_HH
#define DUNE_GRID_TEST_CHECKGEOMETRY_HH

#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/test/checkgeometry.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  template<class Grid>
  struct GeometryChecker;

  /** \param geometry The local geometry to be tested
   * \param type The type of the element that the local geometry is embedded in
   * \param geoName Helper string that will appear in the error message
   */
  template< int mydim, int cdim, class Grid, template< int, int, class > class Imp >
  void checkLocalGeometry ( const Geometry< mydim, cdim, Grid, Imp > &geometry,
                            GeometryType type, const std::string &geoName = "local geometry" )
  {
    checkGeometry( geometry );

    // check that corners are within the reference element of the given type
    assert( type.dim() == cdim );

    // we can't get the reference element from the geometry here, because we need global coordinates
    auto refElement = referenceElement< typename Grid::ctype, cdim >( type );

    const int numCorners = geometry.corners();
    for( int i = 0; i < numCorners; ++i )
    {
      if( !refElement.checkInside( geometry.corner( i ) ) )
      {
        std::cerr << "Corner " << i
                  << " of " << geoName << " not within reference element: "
                  << geometry.corner( i ) << "." << std::endl;
      }
    }
  }


  template <class GI>
  struct CheckSubEntityGeometry
  {
    template <int codim>
    struct Operation
    {
      template <class Entity>
      static void apply(const Entity &entity)
      {
        std::integral_constant<
          bool, Dune::Capabilities::hasEntity<GI,codim>::v
          > capVar;
        check(capVar,entity);
      }
      template<class Geometry>
      static void checkGeometryStatic(const Geometry& geometry)
      {
        std::integral_constant<
          bool, Dune::Capabilities::hasGeometry<GI,codim>::v
          > capVar;
        checkGeometry(capVar, geometry);
      }
      template<class Geometry>
      static void checkGeometry(const std::true_type&,const Geometry& geometry)
      {
        Dune::checkGeometry(geometry);
      }
      template<class Geometry>
      static void checkGeometry(const std::false_type&,const Geometry& /*geometry*/)
      {}
      template <class Entity>
      static void check(const std::true_type&, const Entity &entity)
      {
        for (unsigned int i=0; i<entity.subEntities(codim); ++i)
          {
            const auto subEn = entity.template subEntity<codim>(i);
            auto subGeo = subEn.geometry();

            if( subEn.type() != subGeo.type() )
              std::cerr << "Error: Entity and geometry report different geometry types on codimension " << codim << "." << std::endl;

            // Move from dynamic codim to static codim to
            // prevent checking non-existing geometries
            switch(codim+1)
              {
              case 0:
                {
                  Operation<0> geometryChecker;
                  geometryChecker.checkGeometryStatic(subGeo);
                  break;
                }
              case 1:
                {
                  Operation<1> geometryChecker;
                  geometryChecker.checkGeometryStatic(subGeo);
                  break;
                }
              case 2:
                {
                  Operation<2> geometryChecker;
                  geometryChecker.checkGeometryStatic(subGeo);
                  break;
                }
              case 3:
                {
                  Operation<3> geometryChecker;
                  geometryChecker.checkGeometryStatic(subGeo);
                  break;
                }
              default:
                break;
              }
          }
      }
      template <class Entity>
      static void check(const std::false_type&, const Entity &)
      {}
    };
  };

  namespace Impl
  {

    // Stores copies of geometries together with a snapshot of their corners
    // taken at the moment the geometry was obtained from the grid.
    template< class Geo >
    class GeometryLifetimeSnapshots
    {
      struct Snapshot
      {
        Geo geometry;
        GeometryType type;
        std::vector< typename Geo::GlobalCoordinate > corners;
        std::string description;
      };

    public:
      void push_back ( const Geo &geometry, std::string description )
      {
        Snapshot snapshot{ geometry, geometry.type(), {}, std::move( description ) };
        for( int i = 0; i < geometry.corners(); ++i )
          snapshot.corners.push_back( geometry.corner( i ) );
        snapshots_.push_back( std::move( snapshot ) );
      }

      // check that the stored geometries still describe the snapshot
      void check () const
      {
        using ctype = typename Geo::ctype;
        const ctype tolerance = 8 * std::numeric_limits< ctype >::epsilon();
        for( const Snapshot &snapshot : snapshots_ )
        {
          const Geo &geometry = snapshot.geometry;
          if( (geometry.type() != snapshot.type) || (std::size_t( geometry.corners() ) != snapshot.corners.size()) )
            DUNE_THROW( InvalidStateException, "Stored copy of " << snapshot.description
                        << " changed its type or number of corners after further grid traversal." );
          for( std::size_t i = 0; i < snapshot.corners.size(); ++i )
          {
            const auto corner = geometry.corner( i );
            using std::max;
            if( (corner - snapshot.corners[ i ]).two_norm() > tolerance * max( ctype( 1 ), snapshot.corners[ i ].two_norm() ) )
              DUNE_THROW( InvalidStateException, "Stored copy of " << snapshot.description
                          << " changed after further grid traversal: corner " << i << " is " << corner
                          << " but was " << snapshot.corners[ i ] << "." );
          }
        }
      }

    private:
      std::vector< Snapshot > snapshots_;
    };

    // traverse the grid view and evaluate all geometries to touch any
    // internal caches the grid implementation might have
    template< class GV >
    void touchAllGeometries ( const GV &gridView )
    {
      using Grid = typename GV::Grid;
      for( const auto &element : elements( gridView ) )
      {
        if( element.partitionType() == GhostEntity )
          continue;

        element.geometry().corner( 0 );
        Hybrid::forEach( std::make_index_sequence< GV::dimension >{}, [ & ] ( auto i ) {
          constexpr int codim = i + 1;
          if constexpr (Capabilities::hasEntity< Grid, codim >::v)
            for( unsigned int j = 0; j < element.subEntities( codim ); ++j )
              element.template subEntity< codim >( j ).geometry().corner( 0 );
        } );
        if( element.hasFather() )
          element.geometryInFather().corner( 0 );

        for( const auto &intersection : intersections( gridView, element ) )
        {
          intersection.geometry().corner( 0 );
          intersection.geometryInInside().corner( 0 );
          if( intersection.neighbor() )
            intersection.geometryInOutside().corner( 0 );
        }
      }
    }

  } // namespace Impl


  /** \brief check that geometries remain valid until the grid is modified
   *
   *  The grid interface guarantees that geometry objects returned by
   *  entities and intersections remain valid until the grid is modified (or
   *  deleted). In particular, neither advancing the iterator nor destroying
   *  the entity or intersection the geometry was obtained from may change a
   *  stored copy of the geometry.
   *
   *  This check stores copies of the element geometries, the sub-entity
   *  geometries, the geometries in father, and of all intersection geometries
   *  for the first \p checkElementCount elements. After a further traversal
   *  of the whole grid view, the stored geometries are compared against
   *  their state at the time they were obtained.
   */
  template<typename GV>
  void checkGeometryLifetime (const GV &gridView, std::size_t checkElementCount = 32)
  {
    using Grid = typename GV::Grid;
    using Element = typename GV::template Codim< 0 >::Entity;
    using Intersection = typename GV::Intersection;

    Impl::GeometryLifetimeSnapshots< typename Element::Geometry > elementGeometries;
    Impl::GeometryLifetimeSnapshots< typename Element::LocalGeometry > fatherGeometries;
    Impl::GeometryLifetimeSnapshots< typename Intersection::Geometry > intersectionGeometries;
    Impl::GeometryLifetimeSnapshots< typename Intersection::LocalGeometry > intersectionLocalGeometries;

    auto subEntityGeometryStorage = unpackIntegerSequence( [ & ] ( auto... i ) {
      return std::make_tuple( Impl::GeometryLifetimeSnapshots< typename GV::template Codim< i+1 >::Geometry >{}... );
    }, std::make_index_sequence< GV::dimension >{} );

    std::size_t count = 0;
    for( const auto &element : elements( gridView ) )
    {
      if( count >= checkElementCount )
        break;
      if( element.partitionType() == GhostEntity )
        continue;

      const std::string elementName = "element " + std::to_string( count );
      ++count;

      elementGeometries.push_back( element.geometry(), "geometry of " + elementName );

      Hybrid::forEach( std::make_index_sequence< GV::dimension >{}, [ & ] ( auto i ) {
        constexpr int codim = i + 1;
        if constexpr (Capabilities::hasEntity< Grid, codim >::v)
          for( unsigned int j = 0; j < element.subEntities( codim ); ++j )
            std::get< i >( subEntityGeometryStorage ).push_back( element.template subEntity< codim >( j ).geometry(),
              "geometry of sub-entity " + std::to_string( j ) + " of codimension " + std::to_string( codim ) + " of " + elementName );
      } );

      if( element.hasFather() )
        fatherGeometries.push_back( element.geometryInFather(), "geometryInFather of " + elementName );

      for( const auto &intersection : intersections( gridView, element ) )
      {
        const std::string intersectionName = "intersection " + std::to_string( intersection.indexInInside() ) + " of " + elementName;
        intersectionGeometries.push_back( intersection.geometry(), "geometry of " + intersectionName );
        intersectionLocalGeometries.push_back( intersection.geometryInInside(), "geometryInInside of " + intersectionName );
        if( intersection.neighbor() )
          intersectionLocalGeometries.push_back( intersection.geometryInOutside(), "geometryInOutside of " + intersectionName );
      }
    }

    // traverse the grid again
    Impl::touchAllGeometries( gridView );

    elementGeometries.check();
    Hybrid::forEach( std::make_index_sequence< GV::dimension >{}, [ & ] ( auto i ) {
      std::get< i >( subEntityGeometryStorage ).check();
    } );
    fatherGeometries.check();
    intersectionGeometries.check();
    intersectionLocalGeometries.check();
  }

  template<class Grid>
  struct GeometryChecker
  {
    template<int codim>
    using SubEntityGeometryChecker =
      typename CheckSubEntityGeometry<Grid>::template Operation<codim>;

    template< class VT >
    void checkGeometry ( const GridView< VT > &gridView )
    {
      const auto end = gridView.template end<0>();
      auto it = gridView.template begin<0>();
      for( ; it != end; ++it )
        Hybrid::forEach(std::make_index_sequence<GridView<VT>::dimension+1>{},[&](auto i){SubEntityGeometryChecker<i>::apply(*it);});
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKGEOMETRY_HH
