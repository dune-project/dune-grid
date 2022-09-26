# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import dune.grid

mshfile = "../../../doc/grids/gmsh/circle1storder.msh"
unstructuredGrid = dune.grid.ugGrid( (dune.grid.reader.gmsh, mshfile), dimgrid=2 )
