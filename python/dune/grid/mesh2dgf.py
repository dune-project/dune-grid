# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import numpy as np

def meshDim(cells):
    cell_names = {c.lower() for c in cells}
    three_d = {"tetra", "hexahedron", "pyramid", "wedge", "tetra10", "hexa20"}
    two_d = {"triangle", "quad", "triangle6", "quad9"}
    elementType="general"
    if cell_names & three_d:
        if len(cell_names & three_d) == 1:
            if cell_names & {"tetra"}: elementType="simplex"
            elif cell_names & {"hexahedron"}: elementType="cube"
        return 3, elementType
    if cell_names and (cell_names <= two_d or (cell_names & two_d and not (cell_names & three_d))):
        if len(cell_names & two_d) == 1:
            if cell_names & {"triangle"}: elementType="simplex"
            elif cell_names & {"quad"}: elementType="cube"
        return 2, elementType
    return 1, "simplex"

def mesh2DGF(mesh, bndDomain = None, bndSegments = None, periodic = None, dim = None,
             computeFirstIndex = True, defaultBndId = None, ignoreInternalId=0):
    """
    Parameter:
       mesh        either a filename then meshio is used to read the mesh file
                   or pair of (points,cells) with
                   points      array(list) of points of length dim
                   cells       dict containing element vertex numbers
       bndDomain   dict id -> list[lower,upper] (or id -> str) where lower and upper describe the bounding box of the boundary section
       bndSegments dict id -> list of lists containing vertex numbers of boundary segments
       periodic    string containing periodic boundary transformation
       dim         dimension of grid
       computeFirstIndex if true, first vertex index is computed bases on all cells vertices (default is on)

    Returns
        String containing the mesh description in DGF format.
    """
    if not type(mesh) is str:
        points, cells = mesh
        dim, elementType = meshDim(cells)
    else:
        try:
            import meshio
        except ImportError:
            raise ImportError("Function `mesh2DGF` uses the `meshio` package - run `pip install meshio`")
        mesh = meshio.read(mesh)
        dim, elementType = meshDim(mesh.cells_dict)
        if dim == 2:
            bndCells = ["line"]
        else:
            bndCells = {"triangle", "quad"}

        cells = mesh.cells_dict
        points = mesh.points.astype("float")

        if bndSegments is None:
            bndSegments = {}
        for cell in bndCells:
            # we prefer 'physical' tagging
            try:
                bnd = list(cells[cell])
            except KeyError:
                continue
            try:
                for i,(line,id) in enumerate( zip(bnd,mesh.cell_data_dict["gmsh:physical"][cell]) ):
                    if id == ignoreInternalId: continue # inside skeleton tag
                    if id <= 0: continue
                    if id in bndSegments:
                        bndSegments[id] += [line]
                    else:
                        bndSegments[id] = [line]
                    bnd[i] = None
            except KeyError:
                pass
            # we prefer 'physical' but will try 'geometrical' tagging as well
            try:
                for line,id in zip(bnd,mesh.cell_data_dict["gmsh:geometrical"][cell]):
                    if line is None: continue
                    if id == ignoreInternalId: continue # inside skeleton tag
                    if id <= 0: continue
                    if id in bndSegments:
                        bndSegments[id] += [line]
                    else:
                        bndSegments[id] = [line]
            except KeyError:
                pass

    if defaultBndId:
        assert not bndDomain
        bndDomain = {defaultBndId:"default"}
    if dim is None:
        if "tetra" in cells or "hexahedron" in cells:
            dim = 3
        elif "quad" in cells or "triangle" in cells:
            dim = 2
        else:
            dim = len(points[0])

    simplex = "triangle" if "triangle" in cells else None
    if dim == 3 and "tetra" in cells:
        simplex = "tetra"

    # default numbering is 0 based
    firstVertexIndex = 0
    if computeFirstIndex:
        for key, cellVertices in cells.items():
            firstVertexIndex = min(firstVertexIndex, np.min(cellVertices))

    dgf="DGF\nVertex\n"
    # if first vertex index is not 0 we need to add this here
    if firstVertexIndex > 0:
        dgf += f"firstindex {firstVertexIndex}\n"

    for p in points:
        for i in range(dim):
            dgf += str(p[i]) + " "
        dgf += "\n"
    dgf += "#\n\n"

    if simplex is not None:
        dgf += "Simplex\n"
        for t in cells[simplex]:
            for v in t:
                dgf += str(v) + " "
            dgf += "\n"
        dgf += "#\n\n"

    if "quad" in cells and dim==2:
        dgf += "Cube\n"
        # gmsh has a different reference quadrilateral
        vxmap = [0,1,3,2] # flip vertex 2 and 3
        for t in cells["quad"]:
            for i in range(4):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"
    if "hexahedron" in cells and dim==3:
        dgf += "Cube\n"
        # gmsh has a different reference hexahedron
        vxmap = [0,1,3,2,4,5,7,6] # flip vertex 2,3 and 6,7
        for t in cells["hexahedron"]:
            for i in range(8):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"

    # boundary segments
    if bndSegments is not None:
        assert isinstance(bndSegments, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundarySegments\n"
        for bndid,bndsegs in bndSegments.items():
            for segment in bndsegs:
                dgf += str(bndid)
                for vx in segment:
                    dgf += " " + str(vx)
                dgf += "\n"
        dgf += "#\n\n"

    # boundary domain section
    if bndDomain is not None:
        assert isinstance(bndDomain, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundaryDomain\n"
        for bndid,bnd in bndDomain.items():
            if isinstance(bnd,str):
                if bnd == 'default':
                    dgf += bnd + " " + str(bndid) +"\n"
                else:
                    dgf += str(bndid) + " " + bnd +"\n"
            else: # tuple or list
                # has to be either list or tuple here
                assert isinstance(bnd,(tuple,list))
                assert len(bnd) == 2 # a lower left and upper right corner
                dgf += str(bndid)
                for coord in bnd:
                    assert len(coord) == dim # should a coordinate in the domain
                    for c in coord:
                        dgf += " " + str(c)
                dgf += "\n"

        dgf += "#\n"

    # periodic boundaries
    if periodic is not None:
        dgf += periodic

    return dgf, dim, elementType

def main() -> int:
    import sys
    from dune.alugrid.importmesh import importMesh
    from dune.alugrid import aluGrid as leafGridView
    from dune.fem.function import boundaryFunction
    import matplotlib.pyplot as plt
    from dune.grid import reader
    msh = sys.argv[1]

    if False:
        # test the gmsh reader - boundary ids are missing and 2d cubes fail (3d cubes work...)
        # gridView = leafGridView((reader.gmsh,msh)) # , defaultBndId=7, dimgrid=2, elementType="simplex")
        gridView = leafGridView((reader.gmsh,msh), dimgrid=3, elementType="cube")
        bndIds = boundaryFunction(gridView)
        if gridView.dimension == 2:
            bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
        else:
            gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # plot 1
    gridView = leafGridView((reader.meshio,msh), defaultBndId=5)
    bndIds = boundaryFunction(gridView)
    if gridView.dimension == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # plot 2
    domain = importMesh(msh, defaultBndId=5, ignoreInternalId=5)
    if domain["elementType"] == "general":
        print("can't read in grid with general element type")
        return
    gridView = leafGridView(domain)
    bndIds = boundaryFunction(gridView)
    if gridView.dimension == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # plot 3
    domain, dim, elementType = mesh2DGF(msh, defaultBndId = 5)
    try:
        with open ("test.dgf", "w") as file:
            file.write(domain)
    except:
        pass
    gridView = leafGridView((reader.dgfString,domain), dimgrid=dim, elementType=elementType)
    bndIds = boundaryFunction(gridView)
    if dim == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # plot 4: without a default - this will use '1'
    gridView = leafGridView((reader.meshio,msh))
    bndIds = boundaryFunction(gridView)
    if gridView.dimension == 2:
        bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
    else:
        gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)

    # plot 5
    try: # without a default using DGF - this should give an exception
         # since boundary segs are attached but no default given
        domain, dim, elementType = mesh2DGF(msh)
        gridView = leafGridView((reader.dgfString,domain), dimgrid=dim, elementType=elementType)
        bndIds = boundaryFunction(gridView)
        if gridView.dimension == 2:
            bndIds.plot(level=0, gridLines="white", cmap=plt.cm.jet, block=False)
        else:
            gridView.writeVTK(msh, pointdata=[bndIds], subsampling=0)
    except RuntimeError:
        print("failed due to missing default")

    plt.show()

if __name__ == "__main__":
    raise SystemExit(main())
