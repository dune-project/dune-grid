// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_DATAHANDLE_HH
#define DUNE_GRID_CONCEPTS_DATAHANDLE_HH

#include <concepts>

#include <dune/grid/concepts/archetypes/datahandle.hh>
#include <dune/grid/concepts/archetypes/entity.hh>
#include <dune/grid/concepts/archetypes/messagebuffer.hh>

namespace Dune::Concept {

template <class DH>
concept CommDataHandle = requires(const DH chandle, const Archetypes::Entity<2,0>& entity)
{
  typename DH::DataType;

  { chandle.contains(/*dim*/ 0, /*codim*/ 0)  } -> std::convertible_to<bool>;
  { chandle.fixedSize(/*dim*/ 0, /*codim*/ 0) } -> std::convertible_to<bool>;
  { chandle.size(entity)                      } -> std::integral;

  requires requires(DH handle, Archetypes::MessageBuffer<typename DH::DataType> buffer)
  {
    handle.gather(buffer, entity);
    handle.scatter(buffer, entity, /*size*/ 0u);
  };
};

static_assert(CommDataHandle< Archetypes::CommDataHandle<double> >);

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPTS_DATAHANDLE_HH
