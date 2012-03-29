// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERONED_HH
#define DUNE_DGFPARSERONED_HH

//- C++ includes
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <istream>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-grid includes
#include <dune/grid/common/intersection.hh>
#include <dune/grid/onedgrid.hh>

//- local includes
#include "dgfparser.hh"


namespace
{
  // helper method used below
  double getfirst ( std::vector< double > v )
  {
    assert( v.size() == 1 );
    std::cout << "v[ 0 ] = " << v[ 0 ] << std::endl;
    return v[ 0 ];
  }
}  // end anonymous namespace



namespace Dune
{

  // DGFGridInfo
  // -----------

  template< >
  struct DGFGridInfo< OneDGrid >
  {
    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return 0.5;
    }
  };



  // DGFGridFactory< OneDGrid >
  // --------------------------

  template< >
  struct DGFGridFactory< OneDGrid >
  {
    /** \brief grid type */
    typedef OneDGrid Grid;
    /** \brief grid dimension */
    const static int dimension = Grid::dimension;
    /** \brief MPI communicator type */
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    /** \brief constructor taking istream */
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : grid_( 0 ),
        emptyParameters_( 0 )
    {
      generate( input, comm );
    }

    /** \brief constructor taking filename */
    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : grid_( 0 ),
        emptyParameters_( 0 )
    {
      std::ifstream input( filename.c_str() );
      generate( input, comm );
    }

    /** \brief get grid */
    Grid *grid () const
    {
      return grid_;
    }

    /** \brief always returns false */
    template < class GG, template< class > class II >
    bool wasInserted ( const Dune::Intersection< GG, II > &intersection ) const
    {
      return false;
    }

    template < class GG, template< class > class II >
    int boundaryId ( const Dune::Intersection< GG, II > &intersection ) const
    {
      // OneDGrid returns boundary segment index;
      // we return the index as the method boundaryId is deprecated
      return intersection.boundarySegmentIndex();
    }

    /** \brief OneDGrid does not support parameters, returns 0 */
    template< class Entity >
    int numParameters ( const Entity & ) const
    {
      return 0;
    }

    /** \brief OneDGrid does not support parameters, returns 0 */
    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector< double >& parameter ( const Entity &entity )
    {
      return parameter< Entity::codimension >( entity );
    }

    /** \brief return empty vector */
    template< int codim >
    std::vector< double > &parameter ( const typename Grid::Codim< codim >::Entity &element )
    {
      return emptyParameters_;
    }

    /** \brief OneDGrid does not support boundary parameters */
    bool haveBoundaryParameters () const
    {
      return false;
    }

    /** \brief return invalid default value */
    template < class GG, template< class > class II >
    const DGFBoundaryParameter::type &boundaryParameter ( const Dune::Intersection< GG, II > &intersection ) const
    {
      return DGFBoundaryParameter::defaultValue();
    }

  private:
    // generate grid
    void generate ( std::istream &input, MPICommunicatorType comm );

    Grid *grid_;
    std::vector< double > emptyParameters_;
  };



  // Implementation of DGFGridFactory< OneDGrid >
  // --------------------------------------------

  inline void DGFGridFactory< OneDGrid >::generate ( std::istream &input, MPICommunicatorType comm )
  {
    // try to open interval block
    dgf::IntervalBlock intervalBlock( input );

    // try to open vertex block
    int dimensionworld = Grid::dimensionworld;
    dgf::VertexBlock vertexBlock( input, dimensionworld );

    // check at least one block is active
    if( !( vertexBlock.isactive() || intervalBlock.isactive() ))
      DUNE_THROW( DGFException, "No readable block found" );

    int nvertices = 0;
    std::vector< std::vector< double > > vertices;

    // read vertices first
    if( vertexBlock.isactive() )
    {
      int nparameter = 0;
      std::vector< std::vector< double > > parameter;
      nvertices = vertexBlock.get( vertices, parameter, nparameter );

      if( nparameter > 0 )
        std::cerr << "Warning: vertex parameters will be ignored" << std::endl;
    }

    // get vertices from interval block
    if ( intervalBlock.isactive() )
    {
      if( intervalBlock.dimw() != dimensionworld )
      {
        DUNE_THROW( DGFException, "Error: wrong coordinate dimension in interval block \
                                   (got "                                                                                         << intervalBlock.dimw() << ", expected " << dimensionworld << ")" );
      }

      int nintervals = intervalBlock.numIntervals();
      for( int i = 0; i < nintervals; ++i )
        nvertices += intervalBlock.getVtx( i, vertices );
    }

    // remove duplicates
    std::sort( vertices.begin(), vertices.end() );
    std::vector< std::vector< double > >::iterator it = std::unique( vertices.begin(), vertices.end() );
    vertices.erase( it, vertices.end() );
    if( nvertices != int( vertices.size() ) )
    {
      std::cerr << "Warning: removed duplicate vertices" << std::endl;
      nvertices = vertices.size();
    }

    // copy to vector of doubles
    std::vector< double > vtx( nvertices );
    transform( vertices.begin(), vertices.end(), vtx.begin(), getfirst );

    grid_ = new OneDGrid( vtx );
  }

} // end namespace Dune

#endif // #ifndef DUNE_DGFPARSERONED_HH
