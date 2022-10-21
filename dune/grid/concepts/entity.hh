// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPT_ENTITY_HH
#define DUNE_GRID_CONCEPT_ENTITY_HH


#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/archetypes/entity.hh>

#include <dune/geometry/type.hh>

namespace Dune::Concept {

/**
 * @brief Model of an entity seed
 * @ingroup GridConcepts
 * @details Dune::EntitySeed is a template for this model
 */
template<class S>
concept EntitySeed = requires(S seed)
{
  requires std::default_initializable<S>;
  { S::codimension  } -> std::convertible_to<int>;
  { seed.isValid()  } -> std::convertible_to<bool>;
};

static_assert(EntitySeed< Archetypes::EntitySeed<0> >);


/**
 * @brief Model of a grid entity for any codimension
 * @ingroup GridConcepts
 * @details Dune::Entity is a template for this model.
 */
template<class E>
concept EntityGeneral = requires(E e, unsigned int codim)
{
  requires Geometry<typename E::Geometry>;
  requires EntitySeed<typename E::EntitySeed>;
  requires E::mydimension == (E::dimension - E::codimension);
  { e.level()             } -> std::convertible_to<int>;
  { e.partitionType()     } -> std::convertible_to<Dune::PartitionType>;
  { e.geometry()          } -> std::convertible_to<typename E::Geometry>;
  { e.type()              } -> std::convertible_to<Dune::GeometryType>;
  { e.subEntities(codim)  } -> std::convertible_to<unsigned int>;
  { e.seed()              } -> std::convertible_to<typename E::EntitySeed>;
  { e==e                  } -> std::convertible_to<bool>;
  { e!=e                  } -> std::convertible_to<bool>;
  requires std::default_initializable<E>;
  requires std::copyable<E>;
};

static_assert(EntityGeneral< Archetypes::Entity<2,0> >);


namespace Impl
{
  template<class E, int codim>
  concept EntityCodimExtended = requires(E e, int subEntity)
  {
    { e.template subEntity<codim>(subEntity) } -> EntityGeneral;
  };

  template<class E, std::size_t dim>
  concept AllEntityCodimsExtended = requires(std::make_index_sequence<dim+1> codims)
  {
    []<std::size_t... cc>(std::index_sequence<cc...>)
        requires (Impl::EntityCodimExtended<E,cc> &&...) {} (codims);
  };
}

/**
 * @brief Model of a grid entity with extended requirements for codimension 0
 * @ingroup GridConcepts
 * @details Dune::Entity of codimension 0 is a template for this model.
 */
template<class E>
concept EntityExtended = requires(E e)
{
  requires (E::codimension == 0);
  requires EntityGeneral<E>;
  requires Geometry<typename E::LocalGeometry>;
  { e.father()                      } -> std::convertible_to<E>;
  { e.hasFather()                   } -> std::convertible_to<bool>;
  { e.isLeaf()                      } -> std::convertible_to<bool>;
  { e.isRegular()                   } -> std::convertible_to<bool>;
  { e.geometryInFather()            } -> std::convertible_to<typename E::LocalGeometry>;
  { e.hbegin(/*maxLevel*/ int{})    } -> std::convertible_to<typename E::HierarchicIterator>;
  { e.hend(/*maxLevel*/ int{})      } -> std::convertible_to<typename E::HierarchicIterator>;
  { e.isNew()                       } -> std::convertible_to<bool>;
  { e.mightVanish()                 } -> std::convertible_to<bool>;
  { e.hasBoundaryIntersections()    } -> std::convertible_to<bool>;
  requires Impl::AllEntityCodimsExtended<E, E::dimension>;
  requires std::same_as<E, typename E::template Codim<0>::Entity>;
};

/**
 * @brief Model of a grid entity
 * @ingroup GridConcepts
 * @details Codimension 0 entities are required to have an extended interface.
 * Dune::Entity is a template for this model.
 */
template<class E>
concept Entity = requires {
  requires EntityGeneral<E>;
  requires (E::codimension != 0) || EntityExtended<E>;
};

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPT_ENTITY_HH
