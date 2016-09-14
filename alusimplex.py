from __future__ import absolute_import, division, print_function, unicode_literals

from .create import module

from . import alu

def create(constructor, dimgrid, dimworld=None, refinement="Dune::nonconforming",
        **parameters):

    parameters["type"] = "Dune::simplex"
    parameters["refinement"] = "Dune::nonconforming"
    return alu.create(constructor, dimgrid, dimworld, **parameters)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
