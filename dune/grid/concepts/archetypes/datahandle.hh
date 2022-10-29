// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ARCHETYPES_COMMDATAHANDLE_HH
#define DUNE_GRID_CONCEPTS_ARCHETYPES_COMMDATAHANDLE_HH

#include <cstddef>

#include <dune/grid/common/datahandleif.hh>

#ifndef DOXYGEN
namespace Dune::Concept::Archetypes {

template <class Data>
struct CommDataHandle : public Dune::CommDataHandleIF<CommDataHandle<Data>, Data>
{
  using DataType = Data;

  bool contains (int dim, int codim) const;
  bool fixedSize (int dim, int codim) const;

  template <class Entity>
  std::size_t size (const Entity& entity) const;

  template <class Buffer, class Entity>
  void gather (Buffer& buffer, const Entity& entity) const;

  template <class Buffer, class Entity>
  void scatter (Buffer& buffer, const Entity& entity, std::size_t size);
};

} // end namespace Dune::Concept::Archetypes
#endif // DOXYGEN

#endif // DUNE_GRID_CONCEPTS_ARCHETYPES_COMMDATAHANDLE_HH
