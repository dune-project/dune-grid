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
template<class MB, class T>
concept MessageBuffer = requires(MB mb, T& val)
{
  mb.write(val);
  mb.read(val);
};

static_assert(Concept::MessageBuffer<Archetypes::MessageBuffer<unsigned char>, unsigned char>);

} // namespace Archetypes


} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPT_MESSAGEBUFFER_HH
