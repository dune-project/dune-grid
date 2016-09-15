# from .alberta import create as albertaGrid
# from .alu import create as aluGrid
# from .alusimplex import create as aluSimplexGrid
# from .alucube import create as aluCubeGrid
# from .aluconform import create as aluConformGrid
# from .oned import create as oneDGrid
# from .sp import create as spGrid
# from .ug import create as ugGrid
# from .yasp import create as yaspGrid

from ..common import reader

def cartesianDomain(lower,upper,division,**parameters):
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
    return (reader.dgfString, dgf)

class P1VTKFunction:
    def __init__(self, module, gridView, container):
        self.module = module
        self.mapper = module.MultipleCodimMultipleGeomTypeMapper(gridView, lambda gt: gt.dim == 0)
        self.v = container
    def evaluate(self, e, xi):
        dim = e.dimension
        nVertices = e.subEntities(dim)
        cornerValues = [self.v[self.mapper.subIndex(e, i, dim)] for i in range(nVertices)]
        interpolation = self.module.MultiLinearGeometry(e.type, cornerValues)
        return interpolation.globalPosition(xi)
