// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_DGFEXCEPTION_HH
#define DUNE_GRID_DGFEXCEPTION_HH

#include <dune/common/exceptions.hh>

namespace Dune
{

  //! \brief exception class for IO errors in the DGF parser
  class DGFException
    : public IOError
  {};

}

#endif
