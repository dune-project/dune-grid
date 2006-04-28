// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline YaspGrid<dim,dimworld>*
  MacroGrid ::
  Impl<YaspGrid<dim,dimworld> > ::
  generate(MacroGrid& mg,const char* filename, int MPICOMM) {
    mg.element=Cube;
    std::ifstream gridin(filename);
    IntervalBlock interval(gridin);
    if(!interval.ok()) {
      std::cerr << "Macrofile " << filename << " must have Intervall-Block "
                << "to be used to initialize YaspGrid!\n"
                << "No alternative File-Format defined" << std::endl;
      abort();
    }
    mg.dimw = interval.dimw();
    if (mg.dimw != dimworld) {
      std::cerr << "ERROR: "
                << "Macrofile " << filename << " is for dimension " << mg.dimw
                << " and connot be used to initialize an YaspGrid of dimension "
                << dimworld << std::endl;
      abort();
    }
    FieldVector<double,dimworld> lang;
    FieldVector<int,dimworld>    anz;
    FieldVector<bool,dimworld>   per(false);
    for (int i=0; i<dimworld; i++) {
      lang[i] = interval.length(i);
      anz[i]  = interval.segments(i);
    }
  #if HAVE_MPI_CPP
    if (MPICOMM==-1) {
      std::cerr << "should not be default MPICOMM as last parameter to "
                << "MacroGrid::generate(const char*,int)" << std::endl;
      abort();
    }
    return new YaspGrid<dim,dimworld>(MPICOMM,lang, anz, per , 1 );
  #else
    return new YaspGrid<dim,dimworld>(lang, anz, per , 1 );
  #endif
  }
}
