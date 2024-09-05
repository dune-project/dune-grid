// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPT_ARCHETYPES_ENTITY_HH
#define DUNE_GRID_CONCEPT_ARCHETYPES_ENTITY_HH

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/archetypes/geometry.hh>

#ifndef DOXYGEN
namespace Dune::Concept::Archetypes {

template <int codim>
struct EntitySeed
{
  static constexpr int codimension = codim;
  bool isValid () const;
};


template <int dim, int codim>
struct Entity
{
  static constexpr int dimension = dim;
  static constexpr int codimension = codim;
  static constexpr int mydimension = dim - codim;

  using Geometry = Archetypes::Geometry<mydimension,mydimension>;
  using EntitySeed = Archetypes::EntitySeed<codimension>;

  int level () const;
  Dune::PartitionType partitionType () const;
  Geometry geometry () const;
  Dune::GeometryType type () const;
  unsigned int subEntities (int cd) const;
  EntitySeed seed () const;

  template <int cc>
  Archetypes::Entity<dim,cc> subEntity (int i) const;

  bool operator== (Entity const& entity) const;
  bool operator!= (Entity const& entity) const;
};

} // end namespace Dune::Concept::Archetypes
#endif // DOXYGEN

#endif // DUNE_GRID_CONCEPT_ARCHETYPES_ENTITY_HH
