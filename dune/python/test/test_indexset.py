# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.grid import structuredGrid

def test_subIndices(gridView):
    indexSet = gridView.indexSet
    for intersection in gridView.boundaryIntersections:
        entity      = intersection.inside
        subentity   = (intersection.indexInInside, 1)
        indices_global      = indexSet.subIndices(entity, subentity, 2)
        indices_reference   = entity.referenceElement.subEntities(subentity, 2)
        indices_lookup      = indexSet.subIndices(entity, 2)
        assert len(indices_global) == len(indices_reference)
        for i, j in zip(indices_global, indices_reference):
            assert i == indices_lookup[j]

if __name__ == "__main__":
    gridView = structuredGrid([0,0],[1,1],[10,10])
    test_subIndices(gridView)
