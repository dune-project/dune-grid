// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH

#include <dune/grid/yaspgrid.hh>
#include "dgfparser.hh"

namespace Dune
{

  namespace dgf
  {

    /** \brief Grid parameters for YaspGrid
     *  \ingroup DGFGridParameter
     *
     *  The YaspGridParameter class is in charge of passing parameters specific
     *  to YaspGrid to the grid construction.
     *  Current parameters are:
     *    -# \b overlap defining the overlap of the grid (default value is zero)
     *    .
     *  See the \b examplegrid5.dgf file for an example.
     *
     *  \note The \b periodic parameter has been replaced by the
     *        \b PeriodicFaceTransformation block.
     */
    class YaspGridParameterBlock
      : public GridParameterBlock
    {
    protected:
      int _overlap;     // overlap for YaspGrid

    public:
      //! constructor taking istream
      YaspGridParameterBlock( std::istream &in )
        : GridParameterBlock( in ),
          _overlap( 0 )  // default value
      {
        // check overlap
        if( findtoken( "overlap" ) )
        {
          int x;
          if( getnextentry(x) ) _overlap = x;
          else
          {
            dwarn << "GridParameterBlock: found keyword `overlap' but no value, defaulting to `" <<  _overlap  <<"' !\n";
          }

          if (_overlap < 0)
          {
            DUNE_THROW(DGFException,"Negative overlap specified!");
          }
        }
        else
        {
          dwarn << "YaspGridParameterBlock: Parameter 'overlap' not specified, "
                << "defaulting to '" << _overlap << "'." << std::endl;
        }

      }

      //! get dimension of world found in block
      int overlap () const
      {
        return _overlap;
      }

    };

  }

  template <int dim>
  struct DGFGridFactory< YaspGrid<dim> >
  {
    typedef YaspGrid<dim> Grid;
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

    template <class Intersection>
    bool wasInserted(const Intersection &intersection) const
    {
      return false;
    }
    template <class Intersection>
    int boundaryId(const Intersection &intersection) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector<double> &parameter ( const Entity &entity )
    {
      return emptyParam;
    }

  private:
    void generate( const std::string &filename, MPICommunicatorType comm );
    Grid *grid_;
    std::vector<double> emptyParam;
  };

  template< int dim >
  inline void DGFGridFactory< YaspGrid< dim > >
  ::generate ( const std::string &filename, MPICommunicatorType comm )
  {
    std::ifstream gridin( filename.c_str() );
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

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector<double,dim> lang;
    FieldVector<int,dim>    anz;
    for( int i = 0; i < dim; ++i )
    {
      // check that start point is 0.0
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
      int dir = -1;
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

#if HAVE_MPI
    grid_ = new YaspGrid< dim >( comm, lang, anz, per, grdParam.overlap() );
#else
    grid_ = new YaspGrid< dim >( lang, anz, per, grdParam.overlap() );
#endif
  }

  template <int dim>
  struct DGFGridInfo< YaspGrid<dim> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dim);}
  };


}
#endif
