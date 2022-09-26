// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/vertex.hh>

namespace Dune
{

  namespace dgf
  {

    // VertexBlock
    // -----------

    VertexBlock :: VertexBlock ( std :: istream &in, int &pdimworld )
      : BasicBlock( in, "Vertex" ),
        dimvertex( -1 ),
        dimworld( pdimworld ),
        goodline( true ),
        vtxoffset( 0 ),
        nofParam( 0 )
    {
      if (!isactive())
        return;

      if( findtoken( "firstindex" ) )
      {
        int x;
        if( getnextentry( x ) )
          vtxoffset=x;
      }

      if( findtoken( "parameters" ) )
      {
        int x;
        if( getnextentry( x ) )
          nofParam=x;
      }

      dimvertex = getDimWorld();
      if( pdimworld < 0 )
        pdimworld = dimvertex;
      dimworld = pdimworld;

      if( dimworld < dimvertex )
      {
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": "
                                << "Vertex dimension greater than world dimension." );
      }
      if( dimworld > dimvertex )
      {
        dwarn << id() << " block: Embedding "
              << dimvertex << "-dimensional vertices into "
              << dimworld << "-dimensional space." << std::endl;
      }
    }


    int VertexBlock :: get ( std :: vector< std :: vector< double > > &points,
                             std :: vector< std :: vector< double > > &params,
                             int &nofp )
    {
      nofp = nofParam;
      reset();

      std::vector< double > point( dimworld );
      std::vector< double > param( nofParam );
      while( next( point, param ) )
      {
        points.push_back( point );
        if( nofParam > 0 )
          params.push_back( param );
      }
      return points.size();
    }


    int VertexBlock :: getDimWorld ()
    {
      if( findtoken( "dimension" ) )
      {
        int dimworldFromVertex;
        if( !getnextentry( dimworldFromVertex ) || (dimworldFromVertex <= 0) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid value given for 'dimension'." );
        }
        return dimworldFromVertex;
      }

      reset();
      while( getnextline() )
      {
        int dimworldFromVertex = -nofParam;
        double x;
        while( getnextentry( x ) )
          ++dimworldFromVertex;
        if( dimworldFromVertex > 0 )
          return dimworldFromVertex;
      }

      DUNE_THROW( DGFException,
                  "Error in " << *this << ": "
                              << "Unable to determine dimension of vertices." );
    }


    bool VertexBlock :: next ( std :: vector< double > &point,
                               std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      int n = 0;
      double x;
      for( ; getnextentry( x ); ++n )
      {
        if( n < dimvertex )
          point[ n ] = x;
        else if( n-dimvertex < nofParam )
          param[ n-dimvertex ] = x;
      }

      if( n == 0 )
        return next( point, param );
      else if( n != dimvertex + nofParam )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of coordinates and parameters "
                                               << "(got " << n
                                               << ", expected " << (dimvertex + nofParam) << ")" );
      }

      for( int i = dimvertex; i < dimworld; ++i )
        point[ i ] = double( 0 );
      return (goodline = true);
    }

  } // end namespace dgf

} // end namespace Dune
