// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "nulliteratorfactory.hh"

Dune::OneDGridList<Dune::OneDEntityImp<0> > Dune::OneDGridNullIteratorFactory<0>::emptyList_;
Dune::OneDGridList<Dune::OneDEntityImp<1> > Dune::OneDGridNullIteratorFactory<1>::emptyList_;
