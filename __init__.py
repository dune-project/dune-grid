from .._grid import *
from .core import *

from ._grids import *

from dune.common import FieldVector
from dune.common.compatibility import getNumberOfParameters

registry = dict()

registry["grid"] = grid_registry

def leafGrid(*args, **kwargs):
    return create(*args, **kwargs)

def _getGridFunction(view,y,dimR):
    if dimR is None:
        try:
            dimR = len(y)
        except TypeError:
            dimR = 0
    try:
        GFClass = getattr(view, "GridFunction"+str(dimR))
    except AttributeError:
        raise AttributeError("dune-python not configured to support GridFunctions of dimension "+str(dimR))
    return GFClass

def globalGridFunction(view,GFClass,func):
    def f(element,x):
        return func(element.geometry.position(x))
    def feval(self,element,point=None):
        if point is None:
            return func(element)
        else:
            return func(element.geometry.position(point))
    subclass = type(GFClass.__name__, (GFClass,), {"__call__": feval})
    return subclass(view,f)
def localGridFunction(view,GFClass,func):
    def feval(self,element,point):
        return func(element,point)
    subclass = type(GFClass.__name__, (GFClass,), {"__call__": feval})
    return subclass(view,func)
def gridFunction(view,dimRange=None,isGlobal=None):
    assert hasattr(view, "dimension"), "did you forget to pass in the grid view to the gridFunction decorator"
    def gridFunction_decorator(func):
        if isinstance(isGlobal,bool):
            assert not dimRange is None
            GFClass = _getGridFunction(view,None,dimRange)
            _isGlobal = isGlobal
        else:
            try:
                y = func( FieldVector( [0,]*view.dimension ) )
                _isGlobal = True
            except TypeError:
                try:
                    e = view.elements.__next__()
                except AttributeError:
                    e = view.elements.next()
                y = func( e, FieldVector( [0,]*view.dimension ) )
                _isGlobal = False
            GFClass = _getGridFunction(view,y,dimRange)
        if _isGlobal:
            return globalGridFunction(view,GFClass,func)
        else:
            return localGridFunction(view,GFClass,func)
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
        if not hasattr(cls,"__call__"):
            raise TypeError("Class has no need method")
        if getNumberOfParameters(cls.__call__) == 2: # global case
            return GlobalGridFunction(view)(cls)
        elif getNumberOfParameters(cls.__call__) == 3: # local case
            return LocalGridFunction(view)(cls)
        else:
            raise TypeError("__call__ method needed with 2 or 3 arguments, not %d " %getNumberOfParameters(cls.__call__))
    return GridFunction_decorator
