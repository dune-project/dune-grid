// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH

#include <dune/grid/common/intersection.hh>
#include <dune/grid/yaspgrid.hh>
#include "dgfparser.hh"

namespace Dune
{
  // forward declaration
  // -------------------

  template< class GridImp, template< class > class IntersectionImp >
  class Intersection;


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

  private:
    typedef FieldVector< double, dimension > Point;
    typedef dgf::BoundaryDomBlock BoundaryDomainBlock;

  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      std::ifstream input( filename.c_str() );
      generate( input, comm );
    }

    ~DGFGridFactory ()
    {
      delete boundaryDomainBlock_;
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
      if( boundaryDomainBlock_->isactive() )
      {
        std::vector< Point > corners;
        getCorners( intersection.geometry(), corners );
        const dgf::DomainData *data = boundaryDomainBlock_->contains( corners );
        if( data )
          return data->id();
        else
          return intersection.indexInInside();
      }
      else
        return intersection.indexInInside();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    // return true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return boundaryDomainBlock_->hasParameter();
    }

    template < class GG, template < class > class II >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection< GG, II > & intersection ) const
    {
      if( haveBoundaryParameters() )
      {
        std::vector< Point > corners;
        getCorners( intersection.geometry(), corners );
        const dgf::DomainData *data = boundaryDomainBlock_->contains( corners );
        if( data )
          return data->parameter();
        else
          return DGFBoundaryParameter::defaultValue();
      }
      else
        return DGFBoundaryParameter::defaultValue();
    }

    template< class Entity >
    std::vector<double> &parameter ( const Entity &entity )
    {
      return emptyParam;
    }

  private:
    void generate( std::istream &gridin, MPICommunicatorType comm );

    template< class Geometry >
    static void getCorners ( const Geometry &geometry, std::vector< Point > &corners )
    {
      corners.resize( geometry.corners() );
      for( int i = 0; i < geometry.corners(); ++i )
      {
        const typename Geometry::GlobalCoordinate corner = geometry.corner( i );
        for( int j = 0; j < dimension; ++j )
          corners[ i ][ j ] = corner[ j ];
      }
    }

    Grid *grid_;
    dgf::BoundaryDomBlock *boundaryDomainBlock_;
    std::vector<double> emptyParam;
  };

  template< int dim >
  inline void DGFGridFactory< YaspGrid< dim > >
  ::generate ( std::istream &gridin, MPICommunicatorType comm )
  {
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "YaspGrid can only be created from an interval block." );

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

    boundaryDomainBlock_ = new dgf::BoundaryDomBlock( gridin, dimension );
  }

  template <int dim>
  struct DGFGridInfo< YaspGrid<dim> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dim);}
  };


}
#endif // #ifndef DUNE_DGFPARSERYASP_HH
