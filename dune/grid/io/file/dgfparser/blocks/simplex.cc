// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/simplex.hh>

namespace Dune
{

  namespace dgf
  {

    // SimplexBlock
    // ------------

    SimplexBlock :: SimplexBlock
      ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid )
      : BasicBlock( in, "Simplex" ),
        nofvtx( pnofvtx ),
        vtxoffset( pvtxoffset ),
        dimgrid( pdimgrid ),
        goodline( true ),
        nofparams( 0 )
    {
      if( !isactive() )
        return;

      if( findtoken( "parameters" ) )
      {
        int x = 0;
        if( getnextentry( x ) )
        {
          if( x > 0 )
            nofparams = x;
        }
        if( x <= 0 )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Key 'parameters' found with no or non-positive value." );
        }
      }

      if( dimgrid < 0 )
        dimgrid = getDimGrid();
      pdimgrid = dimgrid;
    }


    int SimplexBlock :: getDimGrid ()
    {
      reset();
      while( getnextline() )
      {
        int count = 0;
        double x;
        while( getnextentry( x ) )
          ++count;
        if( count > nofparams )
          return (count - nofparams) - 1;
      }
      return 0;
    }


    int SimplexBlock
    :: get ( std :: vector< std :: vector< unsigned int > > &simplices,
             std :: vector< std :: vector< double > > &params,
             int &nofp)
    {
      nofp = nofparams;
      reset();

      std :: vector< unsigned int > simplex( dimgrid+1 );
      std :: vector< double > param( nofparams );
      int nofsimpl = 0;
      for( ; next( simplex, param ); ++nofsimpl )
      {
        simplices.push_back( simplex );
        /*
           for( size_t j = 0; j < simplex.size(); ++j )
           simplices[ nofsimpl ][ j ] = simplex[ j ];
         */
        if( nofparams > 0 )
          params.push_back( param );
      }
      return nofsimpl;
    }


    bool SimplexBlock :: next ( std :: vector< unsigned int > &simplex,
                                std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      for( std :: size_t n = 0; n < simplex.size(); ++n )
      {
        int idx;
        if( !getnextentry( idx ) )
        {
          if( n > 0 )
          {
            DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                                   << "Wrong number of vertex indices "
                                                   << "(got " << idx
                                                   << ", expected " << simplex.size() << ")" );
          }
          else
            return next( simplex, param );
        }
        if( (vtxoffset > idx) || (idx >= int(nofvtx + vtxoffset)) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid vertex index "
                                  << "(" << idx << " not in ["
                                  << vtxoffset << ", " << (nofvtx + vtxoffset) << "[)" );
        }
        simplex[ n ] = idx - vtxoffset;
      }

      std :: size_t np = 0;
      double x;
      for( ; getnextentry( x ); ++np )
      {
        if( np < param.size() )
          param[ np ] = x;
      }

      if( np != param.size() )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of simplex parameters "
                                               << "(got " << np
                                               << ", expected " << param.size() << ")" );
      }
      return (goodline = true);
    }


    int SimplexBlock
    :: cube2simplex ( std :: vector< std :: vector< double > > &vtx,
                      std :: vector< std :: vector< unsigned int > > &elements,
                      std :: vector< std :: vector< double > > &params )
    {
      static int offset3[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
                                     {{0,0,0},{1,1,1},{1,0,1},{1,0,0}},
                                     {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
                                     {{0,0,0},{1,1,1},{1,1,0},{0,1,0}},
                                     {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
                                     {{0,0,0},{1,1,1},{0,1,1},{0,0,1}} };
      static int offset2[2][3][2] = {{{0,0},{1,0},{0,1}},
                                     {{1,1},{0,1},{1,0}}};
      if (elements.empty())
        return 0;

      const int dimworld = vtx[ 0 ].size();

      int dimgrid = 0;
      for( size_t n = elements[ 0 ].size(); n > 1; ++dimgrid, n /= 2 ) ;
      if( size_t( 1 << dimgrid ) != elements[ 0 ].size() )
        DUNE_THROW( DGFException, "cube2simplex: all elements must be cubes." );

      dverb << "generating simplices...";
      dverb.flush();

      if( dimgrid == 1 )
        return elements.size();

      std::vector< std::vector< unsigned int > > cubes;
      std::vector< std::vector< double > > cubeparams;
      elements.swap( cubes );
      params.swap( cubeparams );

      if( dimgrid == 3 )
      {
        elements.resize( 6*cubes.size() );
        if( cubeparams.size() > 0 )
          params.resize( 6*cubes.size() );
        for( size_t countsimpl = 0; countsimpl < elements.size(); ++countsimpl )
          elements[ countsimpl ].resize( 4 );
        for( size_t c = 0; c < cubes.size(); ++c )
        {
          for( int tetra = 0; tetra < 6; ++tetra )
          {
            for( int v = 0; v < 4; ++v )
            {
              elements[ c*6+tetra ][ v ]
                = cubes[ c ][ offset3[ tetra ][ v ][ 0 ] +2*offset3[ tetra ][ v ][ 1 ] +4*offset3[ tetra ][ v ][ 2 ] ];
            }
            if( cubeparams.size() > 0 )
              params[ c*6+tetra ] = cubeparams[ c ];
          }
        }
      }
      else if( dimgrid == 2 )
      {
        elements.resize( 2*cubes.size() );
        if( cubeparams.size() > 0 )
          params.resize( 2*cubes.size() );
        for( size_t countsimpl = 0; countsimpl < elements.size(); ++countsimpl )
          elements[ countsimpl ].resize( 3 );
        for( size_t c = 0; c < cubes.size(); ++c )
        {
          int diag = 0;
          double mind = 0;
          for( int d = 0; d < 2; ++d )
          {
            // let's cut the longer diagonal
            double diaglen = 0;
            for( int i = 0; i < dimworld; ++i )
            {
              const double dist = vtx[ cubes[ c ][ d ] ][ i ] - vtx[ cubes[ c ][ 3-d ] ][ i ];
              diaglen += dist*dist;
            }
            if( diaglen < mind )
            {
              mind = diaglen;
              diag = d;
            }
          }
          if( diag == 0 )
          {
            int tmp0 = cubes[ c ][ 0 ];
            cubes[ c ][ 0 ] = cubes[ c ][ 1 ];
            cubes[ c ][ 1 ] = cubes[ c ][ 3 ];
            cubes[ c ][ 3 ] = cubes[ c ][ 2 ];
            cubes[ c ][ 2 ] = tmp0;
          }

          for( int triangle = 0; triangle < 2; ++triangle )
          {
            for( int v = 0; v < 3; ++v )
            {
              elements[ c*2+triangle ][ v ]
                = cubes[ c ][ offset2[ triangle ][ v ][ 0 ] + 2*offset2[ triangle ][ v ][ 1 ] ];
            }
            if( cubeparams.size() > 0 )
              params[ c*2+triangle ] = cubeparams[ c ];
          }
        }
      }
      return elements.size();
    }

  } // end namespace dgf

} // end namespace Dune
