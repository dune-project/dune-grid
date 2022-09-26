// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH

#include <dune/grid/common/intersection.hh>
#include <dune/grid/yaspgrid.hh>
#include "dgfparser.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class GridImp, class IntersectionImp >
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

  /*!
   * \brief Grid factory for YaspGrid with equidistant coordinates.
   */
  template <typename ctype, int dim>
  struct DGFGridFactory< YaspGrid<dim, EquidistantCoordinates<ctype, dim> > >
  {
    typedef YaspGrid<dim, EquidistantCoordinates<ctype, dim> > Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  private:
    typedef FieldVector< ctype, dimension > Point;
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
      if( !input )
        DUNE_THROW( DGFException, "Error: Macrofile '" << filename << "' not found" );
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

    template< class GG, class II >
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
    std::vector<double> &parameter ( const Entity & )
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

  // generate YaspGrid from the provided DGF
  template< typename ctype, int dim >
  inline void DGFGridFactory< YaspGrid< dim , EquidistantCoordinates<ctype, dim> > >
  ::generate ( std::istream &gridin, MPICommunicatorType comm )
  {
    using std::abs;
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "YaspGrid can only be created from an interval block." );

    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "YaspGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dim )
    {
      DUNE_THROW( DGFException,
                  "Cannot read an interval of dimension " << intervalBlock.dimw()
                                                          << " into a YaspGrid< " << dim << " >." );
    }

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector<ctype, dim> lang;
    std::array<int,dim> anz;
    for( int i = 0; i < dim; ++i )
    {
      // check that start point is 0.0
      if( abs( interval.p[ 0 ][ i ] ) > 1e-10 )
      {
        DUNE_THROW( DGFException,
                    "YaspGrid cannot handle grids with non-zero left lower corner." );
      }

      lang[ i ] = interval.p[ 1 ][ i ] - interval.p[ 0 ][ i ];
      anz[ i ]  = interval.n[ i ];
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( gridin, dim );
    std::bitset< dim > per;
    const int numTrafos = trafoBlock.numTransformations();
    for( int k = 0; k < numTrafos; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      bool identity = true;
      for( int i = 0; i < dim; ++i )
        for( int j = 0; j < dim; ++j )
          identity &= (abs( (i == j ? 1.0 : 0.0) - trafo.matrix( i, j ) ) < 1e-10);
      if( !identity )
        DUNE_THROW( DGFException, "YaspGrid can only handle shifts as periodic face transformations." );

      int numDirs = 0;
      int dir = -1;
      for( int i = 0; i < dim; ++i )
      {
        if( abs( trafo.shift[ i ] ) < 1e-10 )
          continue;
        dir = i;
        ++numDirs;
      }
      if( (numDirs != 1) || (abs( abs( trafo.shift[ dir ] ) - lang[ dir ] ) >= 1e-10) )
      {
        std::cerr << "Tranformation '" << trafo
                  << "' does not map boundaries on boundaries." << std::endl;
      }
      else
        per[ dir ] = true;
    }

    // get grid parameters
    dgf::YaspGridParameterBlock grdParam( gridin );

    grid_ = new YaspGrid< dim , EquidistantCoordinates<ctype, dim> >( lang, anz, per, grdParam.overlap(), comm );

    boundaryDomainBlock_ = new dgf::BoundaryDomBlock( gridin, dimension );
  }

  /*!
   * \brief Grid factory for YaspGrid with equidistant coordinates.
   */
  template <typename ctype, int dim>
  struct DGFGridFactory< YaspGrid<dim, EquidistantOffsetCoordinates<ctype, dim> > >
  {
    typedef YaspGrid<dim, EquidistantOffsetCoordinates<ctype, dim> > Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  private:
    typedef FieldVector< ctype, dimension > Point;
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

    template< class GG, class II >
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
    std::vector<double> &parameter ( [[maybe_unused]] const Entity &entity )
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

  // generate YaspGrid from the provided DGF
  template< typename ctype, int dim >
  inline void DGFGridFactory< YaspGrid<dim, EquidistantOffsetCoordinates<ctype, dim> > >
  ::generate ( std::istream &gridin, MPICommunicatorType comm )
  {
    using std::abs;
    dgf::IntervalBlock intervalBlock( gridin );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "YaspGrid can only be created from an interval block." );

    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "YaspGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dim )
    {
      DUNE_THROW( DGFException,
                  "Cannot read an interval of dimension "
                  << intervalBlock.dimw()
                  << " into a YaspGrid< " << dim << " >." );
    }

    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector<ctype, dim> lower;
    FieldVector<ctype, dim> upper;
    std::array<int,dim> anz;
    for( int i = 0; i < dim; ++i )
    {
      lower[ i ] = interval.p[ 0 ][ i ];
      upper[ i ] = interval.p[ 1 ][ i ];
      anz[ i ] = interval.n[ i ];
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( gridin, dim );
    std::bitset< dim > periodic;
    const int numTrafos = trafoBlock.numTransformations();
    for( int k = 0; k < numTrafos; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      bool identity = true;
      for( int i = 0; i < dim; ++i )
        for( int j = 0; j < dim; ++j )
          identity &= (abs( (i == j ? 1.0 : 0.0) - trafo.matrix( i, j ) ) < 1e-10);
      if( !identity )
        DUNE_THROW( DGFException, "YaspGrid can only handle shifts as periodic face transformations." );

      int numDirs = 0;
      int dir = -1;
      for( int currentDir = 0; currentDir < dim; ++currentDir )
      {
        if( abs( trafo.shift[ currentDir ] ) > 1e-10 )
        {
          dir = currentDir;
          ++numDirs;
        }
      }
      if ( (numDirs != 1)
          || (abs( abs( trafo.shift[ dir ] ) - abs( upper[ dir ] - lower[ dir ] ) ) >= 1e-10) )
      {
        std::cerr << "Tranformation '" << trafo
                  << "' does not map boundaries on boundaries." << std::endl;
      }
      else
      {
        periodic[ dir ] = true;
      }
    }

    // get grid parameters
    dgf::YaspGridParameterBlock grdParam( gridin );

    grid_ = new YaspGrid< dim, EquidistantOffsetCoordinates<ctype, dim> >
                        ( lower, upper, anz, periodic, grdParam.overlap(), comm );

    boundaryDomainBlock_ = new dgf::BoundaryDomBlock( gridin, dimension );
  }

  /*!
   * \brief Placeholder for grid factory for YaspGrid with tensor product coordinates.
   *
   * currently tensor product coordinates are currently not supported. Triggers a run time error.
   */
  template< class ctype, int dim >
  class DGFGridFactory< Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> > >
  {
    using Grid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> >;
  public:
    template< typename In >
    DGFGridFactory ( const In & ) {}
    Grid *grid()
    {
      throw std::invalid_argument( "Tensor product coordinates for YaspGrid are currently not supported via the DGFGridFactory" );
    }
  };

  template <typename Coordinates, int dim>
  struct DGFGridInfo< YaspGrid<dim , Coordinates> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return std::pow(0.5,dim);}
  };

}
#endif // #ifndef DUNE_DGFPARSERYASP_HH
