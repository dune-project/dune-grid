# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from ._grid import *
from .core import *

from ._grids import *

from dune.common import FieldVector
from dune.common.utility import getNumberOfParameters

registry = dict()

registry["grid"] = grid_registry

def gridFunction(view,name=None,order=None,dimRange=None):
    assert hasattr(view, "dimension"), "did you forget to pass in the grid view to the gridFunction decorator"
    def gridFunction_decorator(func):
        return view.function(func,name=name,order=order,dimRange=dimRange)
    return gridFunction_decorator
gridFunction._counter = 0

def GridFunction(view, name=None,order=None):
    assert hasattr(view, "dimension"), "did you forget to pass in the grid view to the gridFunction decorator"
    def GridFunction_decorator(cls):
        if not hasattr(cls,"__call__"):
            raise TypeError("Class has no call method")
        class Wrapper(cls):
            def __init__(self, *args, **kwargs):
                cls.__init__(self,*args,**kwargs)
                if getNumberOfParameters(cls.__call__) == 2: # global case
                    self.gf = view.function(lambda x: cls.__call__(self,x),name=name)
                elif getNumberOfParameters(cls.__call__) == 3: # local case
                    self.gf = view.function(lambda e,x: cls.__call__(self,e,x),name=name)
                else:
                    raise TypeError("__call__ method needed with 2 or 3 arguments, not %d " %getNumberOfParameters(cls.__call__))
                if not hasattr(self,"order"):
                    self.order = order
            # note: any magic methods on gf will not be picked up!
            def __getattr__(self, name):
                return getattr(self.gf, name)
            def __call__(self,element,point):
                return self.gf(element,point)
        return Wrapper
    return GridFunction_decorator
