// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{

  template< int dim >
  inline YaspGrid< dim > *
  MacroGrid::Impl< YaspGrid< dim > >
  ::generate ( MacroGrid &mg, const char *filename, MPICommunicatorType MPICOMM )
  {

    std::ifstream gridin( filename );
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
    {
      DUNE_THROW( DGFException,
                  "Macrofile " << filename << " must have Intervall-Block "
                               << "to be used to initialize YaspGrid!\n"
                               << "No alternative File-Format defined");
    }

    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "YaspGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dim )
    {
      DUNE_THROW( DGFException,
                  "Cannot read an interval of dimension " << intervalBlock.dimw()
                                                          << "into a YaspGrid< " << dim << " >." );
    }

    mg.element = Cube;
    mg.dimw = intervalBlock.dimw();

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    // get grid parameters
    dgf::GridParameterBlock grdParam( gridin, true );

    FieldVector<double,dim> lang;
    FieldVector<int,dim>    anz;
    FieldVector<bool,dim>   per(false);

    for( int i = 0; i < dim; ++i )
    {
      // check that start point is > 0.0
      if( fabs( interval.p[ 0 ][ i ] ) > 1e-10 )
      {
        DUNE_THROW( DGFException,
                    "YaspGrid cannot handle grids with non-zero left lower corner." );
      }

      lang[ i ] = interval.p[ 1 ][ i ] - interval.p[ 0 ][ i ];
      anz[ i ]  = interval.n[ i ];
      per[ i ]  = grdParam.isPeriodic( i );
    }

    #if HAVE_MPI
    return new YaspGrid<dim>(MPICOMM,lang, anz, per , grdParam.overlap() );
    #else
    return new YaspGrid<dim>(lang, anz, per , grdParam.overlap() );
    #endif
  }

}
