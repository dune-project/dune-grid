# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

def MultipleCodimMultipleGeomTypeMapper(gridView, layout):
    return gridView.mapper(layout)
