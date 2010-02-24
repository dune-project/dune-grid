// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFS_HH
#define DUNE_DGFS_HH

#include <dune/grid/sgrid.hh>
#include "dgfparser.hh"

namespace Dune
{

  template< int dim, int dimworld, class ctype >
  struct DGFGridFactory< SGrid< dim, dimworld, ctype > >
  {
    typedef SGrid< dim, dimworld, ctype > Grid;

    const static int dimension = Grid::dimension;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      generate( filename, comm );
    }

    Grid *grid() const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      return emptyParam;
    }

  private:
    void generate( const std::string &filename, MPICommunicatorType comm );

    Grid *grid_;
    std::vector< double > emptyParam;
  };



  template< int dim, int dimworld, class ctype >
  inline void DGFGridFactory< SGrid< dim, dimworld, ctype > >
  ::generate ( const std::string &filename, MPICommunicatorType comm )
  {
    std::ifstream gridin( filename.c_str() );
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
    {
      DUNE_THROW( DGFException,
                  "Macrofile " << filename << " must have Intervall-Block "
                               << "to be used to initialize SGrid!\n"
                               << "No alternative File-Format defined" );
    }

    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "SGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dim )
    {
      DUNE_THROW( DGFException,
                  "Cannot read an interval of dimension " << intervalBlock.dimw()
                                                          << "into a SGrid< " << dim << ", " << dimworld << " >." );
    }

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector< double, dimension > lower, upper;
    FieldVector< int, dimension > anz;
    for( int i = 0; i < dimension; ++i )
    {
      lower[ i ] = interval.p[ 0 ][ i ];
      upper[ i ] = interval.p[ 1 ][ i ];
      anz[ i ] = interval.n[ i ];
    }

    grid_ = new Grid( anz, lower, upper );
  }



  template< int dim, int dimworld, class ctype >
  struct DGFGridInfo< SGrid< dim, dimworld, ctype > >
  {
    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return 1.0 / double( 1 << dim );
    }
  };

}
#endif // #ifndef DUNE_DGFS_HH
