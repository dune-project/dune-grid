// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <iostream>

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

    FieldVector<double,dim> lang;
    FieldVector<int,dim>    anz;
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
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( gridin, dim );
    FieldVector< bool, dim > per( false );
    const int numTrafos = trafoBlock.numTransformations();
    for( int k = 0; k < numTrafos; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      bool identity = true;
      for( int i = 0; i < dim; ++i )
        for( int j = 0; j < dim; ++j )
          identity &= (fabs( (i == j ? 1.0 : 0.0) - trafo.matrix( i, j ) ) < 1e-10);
      if( !identity )
        DUNE_THROW( DGFException, "YaspGrid can only handle shifts as periodic face transformations." );

      int numDirs = 0;
      int dir;
      for( int i = 0; i < dim; ++i )
      {
        if( fabs( trafo.shift[ i ] ) < 1e-10 )
          continue;
        dir = i;
        ++numDirs;
      }
      if( (numDirs != 1) || (fabs( fabs( trafo.shift[ dir ] ) - lang[ dir ] ) >= 1e-10) )
      {
        std::cerr << "Tranformation '" << trafo
                  << "' does not map boundaries on boundaries." << std::endl;
      }
      else
        per[ dir ] = true;
    }

    // get grid parameters
    dgf::YaspGridParameterBlock grdParam( gridin );
#if 0
    FieldVector< bool, dim > per( false );
    for( int i = 0; i < dim; ++i )
      per[ i ]  = grdParam.isPeriodic( i );
#endif

    #if HAVE_MPI
    return new YaspGrid<dim>(MPICOMM,lang, anz, per , grdParam.overlap() );
    #else
    return new YaspGrid<dim>(lang, anz, per , grdParam.overlap() );
    #endif
  }

}
