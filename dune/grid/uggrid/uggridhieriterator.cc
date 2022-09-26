// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridhieriterator.hh>

namespace Dune {

//***************************************************************
//
//  --UGGridHierarchicIterator
//  --HierarchicIterator
//
//***************************************************************

template<class GridImp>
void UGGridHierarchicIterator<GridImp>::increment()
{
  if (elementStack_.empty())
    return;

  const int dim = GridImp::dimension;

  const typename UG_NS<dim>::Element* oldTarget = elementStack_.top();
  elementStack_.pop();

  // Traverse the tree no deeper than maxlevel
  if (UG_NS<dim>::myLevel(oldTarget) < maxlevel_) {

    typename UG_NS<dim>::Element* sonList[UG_NS<dim>::MAX_SONS];
    UG_NS<dim>::GetSons(oldTarget, sonList);

    // Load sons of old target onto the iterator stack
    for (int i=0; i<UG_NS<dim>::nSons(oldTarget); i++)
      elementStack_.push(sonList[i]);

  }

  if (elementStack_.empty())
    this->entity_.impl().setToTarget(nullptr,nullptr);
  else
    this->entity_.impl().setToTarget(elementStack_.top(),gridImp_);

}

/////////////////////////////////////////////////////////////////////////////////
//   Explicit template instantiations
/////////////////////////////////////////////////////////////////////////////////

template class UGGridHierarchicIterator<const UGGrid<2> >;
template class UGGridHierarchicIterator<const UGGrid<3> >;

} /* namespace Dune */
