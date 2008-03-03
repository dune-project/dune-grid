// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim>
  inline YaspGrid<dim>*
  MacroGrid ::
  Impl<YaspGrid<dim> > ::
  generate(MacroGrid& mg,const char* filename, MPICommunicatorType MPICOMM)
  {
    mg.element=Cube;
    std::ifstream gridin(filename);
    IntervalBlock interval(gridin);
    if(!interval.isactive()) {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " must have Intervall-Block "
                              << "to be used to initialize YaspGrid!\n"
                              << "No alternative File-Format defined");
    }
    mg.dimw = interval.dimw();

    // get grid parameters
    GridParameterBlock grdParam(gridin, true);

    FieldVector<double,dim> lang;
    FieldVector<int,dim>    anz;
    FieldVector<bool,dim>   per(false);

    for (int i=0; i<dim; i++)
    {
      // check that start point is > 0.0
      if( interval.start(i) < 0.0 )
      {
        DUNE_THROW(DGFException,"YaspGrid cannot handle grids with left lower corner below zero!");
      }

      // set parameter for yaspgrid
      lang[i] = interval.length(i);
      anz[i]  = interval.segments(i);
      per[i]  = grdParam.isPeriodic(i);
    }

  #if HAVE_MPI
    return new YaspGrid<dim>(MPICOMM,lang, anz, per , grdParam.overlap() );
  #else
    return new YaspGrid<dim>(lang, anz, per , grdParam.overlap() );
  #endif
  }
}
