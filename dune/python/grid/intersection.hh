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
        cls.def_property_readonly( "boundary", &Intersection::boundary,
          R"doc(
               Return True if intersection is with interior or exterior boundary.
          )doc"
          );
        cls.def_property_readonly( "boundarySegmentIndex", [] ( const Intersection &i ) {
            return (i.boundary() ? pybind11::cast( i.boundarySegmentIndex() ) : pybind11::none());
          },
          R"doc(
            Index of the boundary segment within the macro grid.

            In many applications, special data needs to be attached to the boundary
            segments of the macro grid (e.g., a function selecting the boundary
            condition).
            Usually, this data is inherited by the children of the boundary segment.

            In the DUNE framework, data is stored in arrays, addressed by an index,
            in this case the boundarySegmentIndex. The size of these arrays can be
            obtained by the Grid::numBoundarySegments method.

            The indices returned by this method are consecutive, zero based, and local to the
            processor. Notice that these indices do not necessarily coincide with the insertion
            indices of the corresponding boundary segments as provided by the GridFactory.
          )doc"
          );

        // conformity
        cls.def_property_readonly( "conforming", &Intersection::conforming,
          R"doc(
            The intersection is called conforming if it is conforming with respect
            to all involved elements, i.e., it does geometrically coincide with the
            indexInInside's facet of the inside element and - if neighbor is true-
            with the indexInOutside's facet of the outside element.

            Returns: True if intersection is conforming.
          )doc"
          );

        // geometric information
        cls.def_property_readonly( "geometry", &Intersection::geometry,
          R"doc(
            Geometrical information about the intersection in global coordinates.

            This method returns a Geometry object that provides a mapping from
            local coordinates of the intersection to global (world) coordinates.

            Note: If the returned geometry has type None then only a limited set of features
                  is available for the geometry, i.e. center and volume.

            Note: The returned geometry object is guaranteed to remain valid until the
                  grid is modified (or deleted).
          )doc"
          );
        cls.def_property_readonly( "type", &Intersection::type,
          R"doc(
            Obtain the type of reference element for this intersection.
          )doc"
          );
        cls.def_property_readonly( "referenceElement", []( const Intersection &self ) {
            return referenceElement< double, mydimension >( self.type() );
          }, pybind11::keep_alive< 0, 1 >(),
          R"doc(
            Obtain the corresponding reference element of this intersection,
            describing the domain of the map.
          )doc" );

        // information on inside entity
        cls.def_property_readonly( "inside", &Intersection::inside,
          R"doc(
            This attribute contains the inside element,
            i.e. the element from which the iteration over intersections was initiated.
          )doc"
          );
        cls.def_property_readonly( "geometryInInside", &Intersection::geometryInInside,
          R"doc(
            Geometrical information about this intersection in local
            coordinates of the inside entity.

            This method returns a Geometry object that provides a mapping from
            local coordinates of the intersection to local coordinates of the
            inside entity.

            Note: The returned geometry object is guaranteed to remain valid until the
                  grid is modified (or deleted).
          )doc"
          );
        cls.def_property_readonly( "indexInInside", &Intersection::indexInInside,
          R"doc(
             Local index of codim 1 entity in the inside entity where
                    intersection is contained in.

             This attribute contains facet number with respect to the generic reference element.

             Note: This index can be used with the inside entity's
                   `subEntity(1)` method to obtain the facet.
          )doc"
          );

        // information on outside entity
        cls.def_property_readonly( "outside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.outside() ) : pybind11::none());
            },
          R"doc(
             This attribute contains the neighboring element if the intersection
             is shared with another element, otherwise this attribute is None .
          )doc"
          );
        cls.def_property_readonly( "geometryInOutside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.geometryInOutside() ) : pybind11::none()); },
          R"doc(
            Geometrical information about this intersection in local
            coordinates of the outside entity.

            This method returns a Geometry object that provides a mapping from
            local coordinates of the intersection to local coordinates of the
            outside entity.

            Note: The returned geometry object is guaranteed to remain valid until the
                  grid is modified (or deleted).

          )doc"
          );
        cls.def_property_readonly( "indexInOutside", [] ( const Intersection &i ) {
              return (i.neighbor() ? pybind11::cast( i.indexInOutside() ) : pybind11::none());
            },
          R"doc(
            Local index of codim 1 entity in outside entity where intersection is contained in.

            This attribute contains the facet number with respect to the generic reference
            element.

            Note: This index can be used with the outside entity's
                  `subEntity(1)` method to obtain the facet.
          )doc"
          );

        // outer normals
        cls.def_property_readonly( "centerUnitOuterNormal", &Intersection::centerUnitOuterNormal,
          R"doc(
            The unit outer normal (length == 1) of the intersection.

            The returned vector is the normal at the center of the
            intersection's geometry. It is scaled to have unit length.
          )doc"
          );
        cls.def( "outerNormal", &Intersection::outerNormal,
          R"doc(
            An outer normal (length not necessarily 1) to the intersection.

            The returned vector may depend on local position within the intersection.
          )doc"
          );
        cls.def( "integrationOuterNormal", &Intersection::integrationOuterNormal,
          R"doc(
            The unit outer normal scaled with the integration element.

            The normal is scaled with the integration element of the intersection. This
            method is redundant but it may be more efficient to use this function
            rather than computing the integration element via geometry().

            The returned vector may depend on local position within the intersection.
          )doc"
          );
        cls.def( "unitOuterNormal", &Intersection::unitOuterNormal,
          R"doc(
            The unit outer normal (length == 1).

            The returned vector may depend on the local position within the intersection.
            It is scaled to have unit length.
          )doc"
          );

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
