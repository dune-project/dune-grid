// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ARCHETYPES_MESSAGEBUFFER_HH
#define DUNE_GRID_CONCEPTS_ARCHETYPES_MESSAGEBUFFER_HH

#ifndef DOXYGEN
namespace Dune::Concept::Archetypes {

template <class DataType>
struct MessageBuffer
{
  void write(const DataType& data);
  void read(DataType& data);
};

} // end namespace Dune::Concept::Archetypes
#endif // DOXYGEN

#endif // DUNE_GRID_CONCEPTS_ARCHETYPES_MESSAGEBUFFER_HH
