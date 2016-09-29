from ._grid import *
from .core import *

from ._grids import *

registry = dict()
registry["grid"] = {
        "Alberta"    : albertaGrid,
        "OneD"       : onedGrid,
        "UG"         : ugGrid,
        "Yasp"       : yaspGrid,
    }

def leafGrid(*args, **kwargs):
    return create(*args, **kwargs)
