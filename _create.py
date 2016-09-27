from ._grids import *

registry = dict()
registry["grid"] = {
        "Alberta"    : albertaGrid,
        "OneD"       : onedGrid,
        "UG"         : ugGrid,
        "Yasp"       : yaspGrid,
    }
