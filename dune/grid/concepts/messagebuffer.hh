// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPT_MESSAGEBUFFER_HH
#define DUNE_GRID_CONCEPT_MESSAGEBUFFER_HH

#include <dune/grid/concepts/archetypes/messagebuffer.hh>

namespace Dune::Concept {

/**
 * @brief Model of a message buffer
 * @ingroup GridConcepts
 */
template<class MB, class DataType>
concept MessageBuffer = requires(MB mb, DataType& data)
{
  mb.write(data);
  mb.read(data);
};

static_assert(Concept::MessageBuffer<Archetypes::MessageBuffer<unsigned char>, unsigned char>);

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPT_MESSAGEBUFFER_HH
