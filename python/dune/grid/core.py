from ._grid import reader
from .map import MultipleCodimMultipleGeomTypeMapper as Mapper

class CartesianDomain(tuple):
    def __new__ (cls, lower,upper,division,**parameters):
        dgf = "DGF\n"
        dgf += "INTERVAL\n"
        dgf += " ".join([str(x) for x in lower]) + "\n"
        dgf += " ".join([str(x) for x in upper]) + "\n"
        dgf += " ".join([str(x) for x in division]) + "\n"
        dgf += "#\n"
        dgf += "GRIDPARAMETER\n"
        for key in parameters:
            dgf += key + " " + str(parameters[key]) + "\n"
        dgf += "#\n"
        try:
            periodic = parameters["periodic"]
            if any(periodic):
                dim = len(lower)
                # create identity matrix in comma separated rows
                eye = ""
                for i in range(dim):
                    for j in range(dim):
                        eye += " 1" if i == j else " 0"
                    eye += " + " if i == dim-1 else ", "

                dgf += "PERIODICFACETRANSFORMATION\n"
                for i,p in enumerate(periodic):
                    if p:
                        dgf += eye
                        dgf += " ".join([str(upper[i]-lower[i]) if i==j
                            else "0" for j in range(dim)])
                        dgf += "\n"
                dgf += "#\n"
        except KeyError:
            pass
        return super(CartesianDomain, cls).__new__(cls,
                       tuple( (reader.dgfString, dgf) ) )
    def __init__(self,lower,upper,division,**parameters):
        self.dimgrid = len(lower)
        self.lower = lower
        self.upper = upper
        self.division = division
        self.param = parameters
    def dimgrid(self):
        return self.dimgrid
def cartesianDomain(lower, upper, division, **parameters):
    return CartesianDomain(lower,upper,division,**parameters)
def structuredGrid(lower,upper,division,**parameters):
    from ._grids import yaspGrid
    domain = cartesianDomain(lower, upper, division, **parameters)
    return yaspGrid(domain, dimgrid=len(lower))

def string2dgf(dgf):
    return (reader.dgfString,"DGF\n" + dgf)


class P1VTKFunction:
    def __init__(self, module, gridView, container):
        self.module = module
        self.mapper = Mapper(gridView, lambda gt: gt.dim == 0)
        self.v = container
    def evaluate(self, e, xi):
        dim = e.dimension
        nVertices = e.subEntities(dim)
        cornerValues = [self.v[self.mapper.subIndex(e, i, dim)] for i in range(nVertices)]
        interpolation = self.module.MultiLinearGeometry(e.type, cornerValues)
        return interpolation.toGlobal(xi)
