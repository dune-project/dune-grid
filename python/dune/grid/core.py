# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from ._grids import CartesianDomain
from .map import MultipleCodimMultipleGeomTypeMapper as Mapper

def cartesianDomain(lower, upper, division, **parameters):
    return CartesianDomain(lower,upper,division,**parameters)
def structuredGrid(lower,upper,division,**parameters):
    from ._grids import yaspGrid
    domain = cartesianDomain(lower, upper, division, **parameters)
    return yaspGrid(domain, dimgrid=len(lower), coordinates="equidistantoffset")

def string2dgf(dgf):
    from ._grid import reader
    return (reader.dgfString,"DGF\n" + dgf)


class P1VTKFunction:
    def __init__(self, module, gridView, container):
        self.module = module
        self.mapper = Mapper(gridView, lambda gt: gt.dim == 0)
        self.v = container
    def evaluate(self, e, xi):
        dim = e.dimension
        nVertices = e.subEntities(dim)
        cornerValues = [self.v[self.mapper.subIndex(e, i, dim)] for i in range(nVertices)]
        interpolation = self.module.MultiLinearGeometry(e.type, cornerValues)
        return interpolation.toGlobal(xi)
