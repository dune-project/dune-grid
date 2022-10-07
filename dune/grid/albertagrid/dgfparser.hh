// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_DGFPARSER_HH
#define DUNE_ALBERTA_DGFPARSER_HH

#include <vector>

#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/blocks/projection.hh>

#include <dune/grid/common/intersection.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class GridImp, class IntersectionImp >
  class Intersection;



  // DGFGridFactory for AlbertaGrid
  // ------------------------------

  template< int dim, int dimworld >
  struct DGFGridFactory< AlbertaGrid< dim, dimworld > >
  {
    typedef AlbertaGrid<dim,dimworld>  Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;
    typedef Dune::GridFactory<Grid> GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );
    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return factory_.wasInserted( intersection );
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return intersection.impl().boundaryId();
    }

    // return true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return dgf_.haveBndParameters;
    }

    template < class GG, class II >
    const DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection< GG, II > & intersection ) const
    {
      typedef Dune::Intersection< GG, II > Intersection;
      typename Intersection::Entity entity = intersection.inside();
      const int face = intersection.indexInInside();

      auto refElem = referenceElement< double, dimension >( entity.type() );
      int corners = refElem.size( face, 1, dimension );
      std :: vector< unsigned int > bound( corners );
      for( int i=0; i < corners; ++i )
      {
        const int k =  refElem.subEntity( face, 1, i, dimension );
        bound[ i ] = factory_.insertionIndex( entity.template subEntity< dimension >( k ) );
      }

      DuneGridFormatParser::facemap_t::key_type key( bound, false );
      const DuneGridFormatParser::facemap_t::const_iterator pos = dgf_.facemap.find( key );
      if( pos != dgf_.facemap.end() )
        return dgf_.facemap.find( key )->second.second;
      else
        return DGFBoundaryParameter::defaultValue();
    }

    template< int codim >
    int numParameters () const
    {
      if( codim == 0 )
        return dgf_.nofelparams;
      else if( codim == dimension )
        return dgf_.nofvtxparams;
      else
        return 0;
    }

    std::vector< double > &parameter ( const Element &element )
    {
      if( numParameters< 0 >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.elParams[ factory_.insertionIndex( element ) ];
    }

    std::vector< double > &parameter ( const Vertex &vertex )
    {
      if( numParameters< dimension >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.vtxParams[ factory_.insertionIndex( vertex ) ];
    }

  private:
    bool generate( std::istream &input );

    Grid *grid_;
    GridFactory factory_;
    DuneGridFormatParser dgf_;
  };



  // DGFGridInfo for AlbertaGrid
  // ---------------------------

  template< int dim, int dimworld >
  struct DGFGridInfo< AlbertaGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
      return dim;
    }

    static double refineWeight ()
    {
      return 0.5;
    }
  };



  // Implementation of DGFGridFactory for AlbertaGrid
  // ------------------------------------------------

  template< int dim, int dimworld >
  inline DGFGridFactory< AlbertaGrid< dim, dimworld > >
  ::DGFGridFactory ( std::istream &input, MPICommunicatorType /* comm */ )
    : dgf_( 0, 1 )
  {
    input.clear();
    input.seekg( 0 );
    if( !input )
      DUNE_THROW(DGFException, "Error resetting input stream." );
    generate( input );
  }


  template< int dim, int dimworld >
  inline DGFGridFactory< AlbertaGrid< dim, dimworld > >
  ::DGFGridFactory ( const std::string &filename, MPICommunicatorType /* comm */ )
    : dgf_( 0, 1 )
  {
    std::ifstream input( filename.c_str() );
    if( !input )
      DUNE_THROW( DGFException, "Macrofile " << filename << " not found." );
    if( !generate( input ) )
      grid_ = new AlbertaGrid< dim, dimworld >( filename.c_str() );
    input.close();
  }

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_DGFPARSER_HH
