from ._grid import *
from .core import *

from ._grids import *

registry = dict()

registry["grid"] = grid_registry

def leafGrid(*args, **kwargs):
    return create(*args, **kwargs)

def globalGridFunction(view):
    def gridFunction_decorator(func):
        def f(element,x):
            return func(element.geometry.position(x))
        gf = view.function(f)
        def feval(self,element,point=None):
            if point:
                return func(element.geometry.position(point))
            else:
                return func(element)
        GFClass = gf.__class__
        subclass = type(GFClass.__name__, (GFClass,), {})
        setattr(subclass, "__call__", feval)
        gf.__class__ = subclass
        return gf
    return gridFunction_decorator
def localGridFunction(view):
    def gridFunction_decorator(func):
        gf = view.function(func)
        def feval(self,element,point):
            return func(element,point)
        GFClass = gf.__class__
        subclass = type(GFClass.__name__, (GFClass,), {})
        setattr(subclass, "__call__", feval)
        gf.__class__ = subclass
        return gf
    return gridFunction_decorator
def gridFunction(view):
    def gridFunction_decorator(func):
        try:
            func( [1./3.,1./3.] )
            return globalGridFunction(view)(func)
        except TypeError:
            return localGridFunction(view)(func)
    return gridFunction_decorator
