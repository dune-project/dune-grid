// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline OneDGrid<dim,dimworld>*
  MacroGrid ::
  Impl<OneDGrid<dim,dimworld> >::generate(MacroGrid& mg,
                                          const char* filename, int MPICOMM) {
    mg.element=Cube;
    std::ifstream gridin(filename);
    if (mg.readDuneGrid(gridin) == 1) {
      std::vector<double> vtxlist(mg.vtx.size());
      for (int i=0; i<vtxlist.size(); ++i)
        vtxlist[i] = mg.vtx[i][0];
      return new OneDGrid<dim,dimworld>(vtxlist);
    }
  }
}
