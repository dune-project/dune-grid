from __future__ import absolute_import, division, print_function, unicode_literals

from dune.common.checkconfiguration import have

def onedGrid(constructor):
    from .grid_generator import module, getDimgrid

    typeName = "Dune::OneDGrid"
    includes = ["dune/grid/onedgrid.hh", "dune/grid/io/file/dgfparser/dgfoned.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

def yaspGrid(constructor, dimgrid=None, coordinates="Dune::EquidistantCoordinates", ctype="double"):
    from .grid_generator import module, getDimgrid

    if not dimgrid:
        dimgrid = getDimgrid(constructor)
    if coordinates == "equidistant":
        coordinates = "Dune::EquidistantCoordinates"
    elif coordinates == "tensorproduct":
        coordinates = "Dune::TensorProductCoordinates"

    typeName = "Dune::YaspGrid< " + str(dimgrid) + ", " + coordinates + "< " + ctype + ", " + str(dimgrid) + " > >"
    includes = ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

grid_registry = {
        "OneD"       : onedGrid,
        "Yasp"       : yaspGrid,
    }

if have("ALBERTA",fail=False):
  def albertaGrid(constructor, dimgrid=None):
      from .grid_generator import module, getDimgrid

      if not dimgrid:
          dimgrid = getDimgrid(constructor)
      typeName = "Dune::AlbertaGrid< " + str(dimgrid) + " >"
      includes = ["dune/grid/albertagrid.hh", "dune/grid/albertagrid/dgfparser.hh"]
      gridModule = module(includes, typeName)

      return gridModule.reader(constructor).leafView

  grid_registry["Alberta"] = albertaGrid

if have("HAVE_DUNE_UGGRID", fail=False):
  def ugGrid(constructor, dimgrid=None, **parameters):
      from .grid_generator import module, getDimgrid

      if not dimgrid:
          dimgrid = getDimgrid(constructor)
      typeName = "Dune::UGGrid< " + str(dimgrid) + " >"
      includes = ["dune/grid/uggrid.hh", "dune/grid/io/file/dgfparser/dgfug.hh"]
      gridModule = module(includes, typeName)

      return gridModule.reader(constructor).leafView

  grid_registry["UG"] = ugGrid

if __name__ == "__main__":
    import doctest
