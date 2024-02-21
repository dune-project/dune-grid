// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_HH

#include <concepts>
#include <cstddef>

#include <dune/geometry/type.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>

namespace Dune::Concept {

/**
 * @brief Model of an intersection
 * @ingroup GridConcepts
 * @details Dune::Grid::Intersection is a template for this model
 */
template<class I>
concept Intersection = std::regular<I> &&
  EntityGeneral<typename I::Entity> &&
  Geometry<typename I::Geometry> &&
  Geometry<typename I::LocalGeometry> &&
requires(const I i, typename I::LocalCoordinate local)
{
  typename I::ctype;
  { I::mydimension                  } -> std::convertible_to<int>;
  { I::dimensionworld               } -> std::convertible_to<int>;
  { i.boundary()                    } -> std::convertible_to<bool>;
  { i.boundarySegmentIndex()        } -> std::convertible_to<std::size_t>;
  { i.neighbor()                    } -> std::convertible_to<bool>;
  { i.inside()                      } -> std::same_as<typename I::Entity>;
  { i.outside()                     } -> std::same_as<typename I::Entity>;
  { i.conforming()                  } -> std::convertible_to<bool>;
  { i.geometryInInside()            } -> std::same_as<typename I::LocalGeometry>;
  { i.geometryInOutside()           } -> std::same_as<typename I::LocalGeometry>;
  { i.geometry()                    } -> std::same_as<typename I::Geometry>;
  { i.type()                        } -> std::same_as<Dune::GeometryType>;
  { i.indexInInside()               } -> std::convertible_to<int>;
  { i.indexInOutside()              } -> std::convertible_to<int>;
  { i.outerNormal(local)            } -> std::convertible_to<typename I::GlobalCoordinate>;
  { i.integrationOuterNormal(local) } -> std::convertible_to<typename I::GlobalCoordinate>;
  { i.unitOuterNormal(local)        } -> std::convertible_to<typename I::GlobalCoordinate>;
  { i.centerUnitOuterNormal()       } -> std::convertible_to<typename I::GlobalCoordinate>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INTERSECTION_HH
