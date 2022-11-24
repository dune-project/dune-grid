# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#######################################################################################################
# This paraview reader adds support for 'dune binary format' files (dbf).
# The file is assumed to be written using 'dune.common.pickle.dump'.
# It therefore consists of two parts (the required jit module source code and
# a pickled list of objects). The list is searched for objects containing
# a 'gridView' attribute - these are all assumed to be grid functions
# over the same grid view and with a 'pointData' attribute.
# If no entry in the list with a 'gridView' attribute is found the first
# entry is assumed to be a grid view and only the grid is plotted.
#
# Additional features:
# --------------------
# - Only pointdata is extracted at the moment - vector data or cell data # would be nice.
# - We always use the dune subsampler so that a non-connected simplex grid is produced.
#   An option of working on the actual grid would be nice to have
# - Add a indicator based grid refinement (or refinement is some part of the domain only)
# - General file reader, i.e., dgf or using GMshReader - issue is figuring out the right dim/dimw
#######################################################################################################

import numpy as np
import os,sys,vtk
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

# In older paraview versions there is no way to set the
# virtual environment to use - use a environment variable
# to set it before starting paraview.
# This should be improved to take other usages into account.

# This finds all 'egg-link' files in a given folder structure.
# These correspond to packages installed 'editable' and need to
# be added by hand to the Python search path:
def find_egglinks(directory_name):
    dune_found = []
    for path, subdirs, files in os.walk(directory_name):
        if not path.endswith("site-packages"):
            continue
        dune_found.append(path)
        for name in files:
            if not "dune" in name:
                continue
            ext = os.path.splitext(name)[1]
            if ext == ".egg-link":
                file_path = os.path.join(path,name)
                with open(file_path,"r") as f:
                    dune_found.append(f.read().split()[0])
    return dune_found
# We find an active virtual env by checking if the environment variable
# 'VIRTUAL_ENV' is set - this work at least with activated environments
# setup with 'venv' on Linux:
def setDuneModulePaths():
    try:
        envdir = os.path.realpath(os.environ['VIRTUAL_ENV'])
        dunePaths = find_egglinks(os.path.join(envdir,"lib"))
        sys.path += dunePaths
        if not "DUNE_PY_DIR" in os.environ:
            os.environ["DUNE_PY_DIR"] = os.path.join(envdir,".cache")
        sys.path += os.path.join(os.environ["DUNE_PY_DIR"],"python","dune","generated")
        # print(os.environ["DUNE_PY_DIR"], dunePaths)
    except KeyError:
        # print("no virtual env path found!")
        pass

##############################################
# Actual reader

dune_extensions = ["dbf"] # ,"dgf"]
@smproxy.reader(
    label="Dune Reader",
    extensions=dune_extensions,
    file_description="dune-supported files",
)
class DuneReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._level = 0
        self._gridView = None
        setDuneModulePaths()
        # do not yet know how to fix the problem in paraview that
        # dune.generated in not found - the following hack works for some reason
        import dune.common
        dune.common.FieldVector([1])
        import dune.common.pickle
        self.load = dune.common.pickle.load

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=dune_extensions, file_description="dune supported files"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.intvector(name="Level", default_values="0")
    @smdomain.intrange(min=0, max=5)
    def SetLevel(self, level):
        self._level = level
        self.Modified()

    def loadData(self):
        ext = os.path.splitext(self._filename)[1]
        if ext == ".dgf":
            print("Still need to implement dgf reading")
            print("Which grid to use with which dimensions?")
        else:
            with open(self._filename,"rb") as f:
                df = self.load(f)
            self._df = [d for d in df if hasattr(d,"gridView")]
            if len(self._df) > 0:
                self._gridView = self._df[0].gridView
            else:
                self._gridView = df[0]
            # make some checks:
            assert hasattr(self._gridView,"dimension"), "file read contains no valid grid view"
            assert all( [hasattr(d,"pointData") for d in self._df] ), "found a non valid grid function (no 'pointData' attribute"
            assert all( [self._gridView == d.gridView for d in self._df] ), "all grid function must be over the same gridView"

    def RequestData(self, request, inInfo, outInfo):
        if self._gridView is None:
            self.loadData()

        points, cells = self._gridView.tessellate(self._level)
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))

        # points need to be 3d:
        if self._gridView.dimWorld == 2:
            vtk_type = vtk.VTK_TRIANGLE
            points = np.hstack([points, np.zeros((len(points), 1))])
        elif self._gridView.dimWorld == 3:
            if self._gridView.dimGrid == 2:
                vtk_type = vtk.VTK_TRIANGLE
            else:
                vtk_type = vtk.VTK_TETRA
        output.SetPoints(points)

        cell_types = np.array([], dtype=np.ubyte)
        cell_offsets = np.array([], dtype=int)
        cell_conn = np.array([], dtype=int)
        ncells, npoints = cells.shape
        cell_types = np.hstack(
                       [cell_types, np.full(ncells, vtk_type, dtype=np.ubyte)]
                     )
        offsets = len(cell_conn) + (1 + npoints) * np.arange(ncells, dtype=int)
        cell_offsets = np.hstack([cell_offsets, offsets])
        conn = np.hstack(
                   [npoints * np.ones((ncells, 1), dtype=int), cells]
               ).flatten()
        cell_conn = np.hstack([cell_conn, conn])
        output.SetCells(cell_types, cell_offsets, cell_conn)  # cell connectivities

        # data
        for df in self._df:
            array = df.pointData(self._level)
            output.PointData.append(array, df.name)  # point data
            # output.CellData.append(array, name)  # cell data
            # output.FieldData.append(array, name)  # field data

        return 1
