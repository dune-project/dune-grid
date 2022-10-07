// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_INTERSECTION_HH
#define DUNE_PYTHON_GRID_INTERSECTION_HH

#include <dune/common/visibility.hh>

#include <dune/python/grid/geometry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/extensions.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // registerGridIntersection
      // ------------------------

      template< class Intersection >
      void registerGridIntersection ( pybind11::handle scope,
          pybind11::class_<Intersection> cls )
      {
        const int mydimension = Intersection::mydimension;
        // information on boundary intersections
        cls.def_property_readonly( "boundary", &Intersection::boundary );
        cls.def_property_readonly( "boundarySegmentIndex", [] ( const Intersection &i ) {
            return (i.boundary() ? pybind11::cast( i.boundarySegmentIndex() ) : pybind11::none());
          } );

        // conformity
        cls.def_property_readonly( "conforming", &Intersection::conforming );

        // geometric information
        cls.def_property_readonly( "geometry", &Intersection::geometry );
        cls.def_property_readonly( "type", &Intersection::type );
        cls.def_property_readonly( "referenceElement", []( const Intersection &self ) {
            return referenceElement< double, mydimension >( self.type() );
          }, pybind11::keep_alive< 0, 1 >(),
          R"doc(
            corresponding reference element, describing the domain of the map
          )doc" );

        // information on inside entity
        cls.def_property_readonly( "inside", &Intersection::inside );
        cls.def_property_readonly( "geometryInInside", &Intersection::geometryInInside );
        cls.def_property_readonly( "indexInInside", &Intersection::indexInInside );

        // information on ouside entity
        cls.def_property_readonly( "outside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.outside() ) : pybind11::none());
            } );
        cls.def_property_readonly( "geometryInOutside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.geometryInOutside() ) : pybind11::none());
            } );
        cls.def_property_readonly( "indexInOutside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.indexInOutside() ) : pybind11::none());
            } );

        // outer normals
        cls.def_property_readonly( "centerUnitOuterNormal", &Intersection::centerUnitOuterNormal);
        cls.def( "outerNormal", &Intersection::outerNormal );
        cls.def( "integrationOuterNormal", &Intersection::integrationOuterNormal );
        cls.def( "unitOuterNormal", &Intersection::unitOuterNormal );

        // comparison
        cls.def( pybind11::self == pybind11::self );
        cls.def( pybind11::self != pybind11::self );
      }

    } // namespace detail



    // registerGridIntersection
    // ------------------------

    template< class GridView >
    auto registerGridIntersection ( pybind11::handle scope )
    {
      typedef typename GridView::Intersection Intersection;
      pybind11::class_< Intersection > cls = insertClass<Intersection>(scope, "Intersection",
         GenerateTypeName(scope,"Intersection")).first;
      detail::registerGridIntersection( scope, cls );
      registerGridGeometry<Intersection>( cls );
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif //#ifndef DUNE_PYTHON_GRID_INTERSECTION_HH
