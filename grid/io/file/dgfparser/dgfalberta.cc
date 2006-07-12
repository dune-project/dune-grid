// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline AlbertaGrid<dim,dimworld>*
  MacroGrid::Impl<AlbertaGrid<dim,dimworld> >::generate(MacroGrid& mg,
                                                        const char* filename,MPI_Comm) {
    mg.element=Simplex;
    std::ifstream gridin(filename);

    std::string str(filename);

    if (mg.readDuneGrid(gridin) == 1) {
      if (mg.dimw != dimworld) {
        std::cerr << "ERROR: "
                  << "Macrofile " << filename << " is for dimension " << mg.dimw
                  << " and connot be used to initialize an "
                  << "AlbertaGrid of dimension "
                  << dimworld << std::endl;
        abort();
      }
      if (mg.dimw == 2) {
        mg.setRefinement(mg.dimw);
        mg.setOrientation(mg.dimw);
      } else {}
      str+=".albertagrid";
      std::ofstream out(str.c_str());
      mg.writeAlberta(out);
    }

    return new AlbertaGrid<dim,dimworld>(str.c_str());
  }
}
