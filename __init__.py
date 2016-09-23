from .grid_generator import create
from .core import *
from .. import alugrid

def leafGrid(*args, **kwargs):
    return create(*args, **kwargs)
