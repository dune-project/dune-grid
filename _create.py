from .alberta import create as albertaGrid
from .oned import create as onedGrid
from .sp import create as spGrid
from .ug import create as ugGrid
from .yasp import create as yaspGrid

registry = dict()
registry["grid"] = {
        "Alberta"    : albertaGrid,
        "OneD"       : onedGrid,
        "SP"         : spGrid,
        "UG"         : ugGrid,
        "Yasp"       : yaspGrid,
    }
