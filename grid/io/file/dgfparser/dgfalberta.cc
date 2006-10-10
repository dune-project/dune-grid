// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline AlbertaGrid<dim,dimworld>*
  MacroGrid::Impl<AlbertaGrid<dim,dimworld> >::generate(MacroGrid& mg,
                                                        const char* filename, MPICommunicatorType ) {
    mg.element=Simplex;
    std::ifstream gridin(filename);

    std::string str(filename);

    if(mg.readDuneGrid(gridin))
    {
      if (mg.dimw != dimworld) {
        DUNE_THROW(DGFException,
                   "Macrofile " << filename << " is for dimension " << mg.dimw
                                << " and connot be used to initialize an "
                                << "AlbertaGrid of dimension "
                                << dimworld);
      }
      mg.setOrientation(0,1);
      mg.setRefinement(0,1,-1,-1);
      str+=".albertagrid";
      std::ofstream out(str.c_str());
      mg.writeAlberta(out);
    }
    return new AlbertaGrid<dim,dimworld>(str.c_str());
  }
}
