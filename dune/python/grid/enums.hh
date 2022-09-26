// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_ENUMS_HH
#define DUNE_PYTHON_GRID_ENUMS_HH

namespace Dune
{

  namespace Python
  {

    enum class Reader { dgf, dgfString, gmsh, structured };
    enum class VTKDataType { CellData, PointData, CellVector, PointVector };
    enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_ENUMS_HH
