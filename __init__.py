from ._grid import *
from .core import *

from ._grids import *

from dune.common import FieldVector
from dune.common.compatibility import getNumberOfParameters

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
            if point is None:
                return func(element)
            else:
                return func(element.geometry.position(point))
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
            func( FieldVector( [0,]*view.dimension ) )
            return globalGridFunction(view)(func)
        except TypeError:
            return localGridFunction(view)(func)
    return gridFunction_decorator

def LocalGridFunction(view):
    def LocalGridFunction_decorator(cls):
        class Wrapper(cls):
            def __init__(self, *args, **kwargs):
                cls.__init__(self,*args,**kwargs)
                self.gf = view.function(self)
            # note: any magic methods on gf will not be picked up!
            def __getattr__(self, name):
                return getattr(self.gf, name)
        return Wrapper
    return LocalGridFunction_decorator
def GlobalGridFunction(view):
    def GlobalGridFunction_decorator(cls):
        class Wrapper(cls):
            def __init__(self, *args, **kwargs):
                cls.__init__(self,*args,**kwargs)
                self.gf = view.function(self)
            # note: any magic methods on gf will not be picked up!
            def __getattr__(self, name):
                return getattr(self.gf, name)
            def __call__(self,element,point=None):
                if point is None:
                    return cls.__call__(self,element)
                else:
                    return cls.__call__(self,element.geometry.position(point))
        return Wrapper
    return GlobalGridFunction_decorator
def GridFunction(view):
    def GridFunction_decorator(cls):
        if getNumberOfParameters(cls.__call__) == 2: # global case
            return GlobalGridFunction(view)(cls)
        elif getNumberOfParameters(cls.__call__) == 3: # local case
            return LocalGridFunction(view)(cls)
        else:
            TypeError("Can't determin correct decorator")
    return GridFunction_decorator
