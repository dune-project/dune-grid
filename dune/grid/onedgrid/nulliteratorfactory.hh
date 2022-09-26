// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_NULL_ITERATORS_HH
#define DUNE_ONEDGRID_NULL_ITERATORS_HH

#include "onedgridlist.hh"

namespace Dune {

  template <int mydim> class OneDEntityImp;

  template <int dim>
  class OneDGridNullIteratorFactory {};

  template <>
  class OneDGridNullIteratorFactory<0> {

  public:

    static OneDGridList<OneDEntityImp<0> >::iterator null() {
      return emptyList_.end();
    }

  private:
    static OneDGridList<OneDEntityImp<0> > emptyList_;
  };

  template <>
  class OneDGridNullIteratorFactory<1> {

  public:

    static OneDGridList<OneDEntityImp<1> >::iterator null() {
      return emptyList_.end();
    }

  private:
    static OneDGridList<OneDEntityImp<1> > emptyList_;
  };

}

#endif
