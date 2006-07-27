// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline SGrid<dim,dimworld>*
  MacroGrid ::
  Impl< SGrid<dim,dimworld> > ::
  generate(MacroGrid& mg,const char* filename, MPICommunicatorType ) {
    mg.element=Cube;
    std::ifstream gridin(filename);
    IntervalBlock interval(gridin);
    if(!interval.ok()) {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " must have Intervall-Block "
                              << "to be used to initialize SGrid!\n"
                              << "No alternative File-Format defined");
    }
    mg.dimw = interval.dimw();
    if (mg.dimw != dimworld) {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " is for dimension " << mg.dimw
                              << " and connot be used to initialize an SGrid of dimension "
                              << dimworld);
    }
    FieldVector<double,dimworld> start;
    FieldVector<double,dimworld> lang;
    FieldVector<int,dimworld>    anz;
    for (int i=0; i<dimworld; i++) {
      start[i] = interval.start(i);
      lang[i] = interval.length(i);
      anz[i]  = interval.segments(i);
    }
    return new SGrid<dim,dimworld>(anz,start,lang);
  }
}
