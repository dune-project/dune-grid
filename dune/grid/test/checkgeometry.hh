// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKGEOMETRY_HH
#define DUNE_GRID_TEST_CHECKGEOMETRY_HH

#include <limits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/test/checkgeometry.hh>

#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridview.hh>

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
      static void checkGeometry(const std::false_type&,const Geometry& geometry)
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

  template<typename GV>
  void checkGeometryLifetime (const GV &gridView)
  {
    typedef typename GV::ctype ctype;
    constexpr static int dim  = GV::dimension;

    const FieldVector<ctype, dim> pos(0.2);

    auto it = gridView.template begin<0>();
    const auto end = gridView.template end<0>();

    // check that it != end otherwise the following is not valid
    if( it == end ) return ;

    #ifndef NDEBUG
    constexpr static int dimw = GV::dimensionworld;
    const FieldVector<ctype, dimw> glob = it->geometry().global(pos);
    #endif

    const auto geomCopy = it->geometry();
    checkGeometry ( geomCopy );

    for( ; it != end; ++it )
    {
      // due to register/memory differences we might have
      // errors < 1e-16
      assert (std::abs((geomCopy.global(pos) - glob).one_norm()) < std::numeric_limits<ctype>::epsilon());
    }
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
