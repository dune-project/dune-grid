// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{

  template< int dim, int dimworld >
  inline SGrid< dim, dimworld > *
  MacroGrid::Impl< SGrid< dim, dimworld > >
  ::generate ( MacroGrid &mg, const char *filename, MPICommunicatorType )
  {
    std::ifstream gridin( filename );
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
    {
      DUNE_THROW( DGFException,
                  "Macrofile " << filename << " must have Intervall-Block "
                               << "to be used to initialize SGrid!\n"
                               << "No alternative File-Format defined");
    }

    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "SGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dimworld )
    {
      DUNE_THROW( DGFException,
                  "Cannot read an interval of dimension " << intervalBlock.dimw()
                                                          << "into a SGrid< " << dim << ", " << dimworld << " >." );
    }

    mg.element = Cube;
    mg.dimw = intervalBlock.dimw();

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector< double, dimworld > lower, upper;
    FieldVector< int, dimworld > anz;

    for( int i = 0; i < dimworld; ++i )
    {
      lower[ i ] = interval.p[ 0 ][ i ];
      upper[ i ] = interval.p[ 1 ][ i ];
      anz[ i ] = interval.n[ i ];
    }

    // SGrid gets the following arguments:
    // - number of cells in each dircetion
    // - position of the lower left corner of the cube
    // - position of the upper right corner of the cube
    return new SGrid< dim, dimworld >( anz, lower, upper );
  }

}
