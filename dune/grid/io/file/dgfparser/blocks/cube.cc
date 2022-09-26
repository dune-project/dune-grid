// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/cube.hh>

namespace Dune
{

  namespace dgf
  {

    // CubeBlock
    // ---------

    CubeBlock :: CubeBlock
      ( std::istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid )
      : BasicBlock( in, "Cube" ),
        nofvtx(pnofvtx),
        dimgrid( pdimgrid ),
        goodline(true),
        map(0),
        nofparams(0),
        vtxoffset(pvtxoffset)
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

      map.resize( 1 << dimgrid );
      for( size_t i = 0; i < map.size(); ++i )
        map[ i ] = i;

      if( findtoken( "map" ) )
      {
        for( size_t i = 0; i < map.size(); ++i )
        {
          int x;
          if( !getnextentry( x ) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": "
                                    << "Incomplete reference mapping "
                                    << "(got " << i << " entries, "
                                    << "expected " << map.size() << " entries." );
          }
          map[ i ] = x;
        }
      }
    }


    int CubeBlock :: getDimGrid ()
    {
      reset();
      while( getnextline() )
      {
        int count = 0;
        double x;
        while( getnextentry( x ) )
          ++count;
        if( count > nofparams )
        {
          count -= nofparams;
          // int dim = (int)(log( count ) / M_LN2);
          int dim = 1;
          while (1<<dim < count)
            dim++;
          if( (dim < 0) || ((1 << dim) != count) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": Number of vertex indices ("
                                    << count << ") is not a power of 2." );
          }
          return dim;
        }
      }
      return 0;
    }


    int CubeBlock :: get ( std :: vector< std :: vector< unsigned int> > &cubes,
                           std :: vector< std :: vector< double > > &params,
                           int &nofp )
    {
      nofp = nofparams;
      reset();

      std :: vector< unsigned int > cube( 1 << dimgrid );
      std :: vector< double > param( nofparams );
      int nofcubes = 0;
      for( ; next( cube, param ); ++nofcubes )
      {
        cubes.push_back( cube );
        /*
           for( size_t j = 0; j < cube.size(); ++j )
           cubes[ nofcubes ][ j ] = cubes[ j ];
         */
        if( nofparams > 0 )
          params.push_back( param );
      }
      return nofcubes;
    }


    bool CubeBlock :: next ( std :: vector< unsigned int > &cube,
                             std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      for( std :: size_t n = 0; n < cube.size(); ++n )
      {
        int idx;
        if( !getnextentry( idx ) )
        {
          if( n > 0 )
          {
            DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                                   << "Wrong number of vertex indices "
                                                   << "(got " << idx
                                                   << ", expected " << cube.size() << ")" );
          }
          else
            return next( cube, param );
        }
        if( (vtxoffset > idx) || (idx >= int(nofvtx + vtxoffset)) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid vertex index "
                                  << "(" << idx << " not in ["
                                  << vtxoffset << ", " << (nofvtx + vtxoffset) << "[)" );
        }
        cube[ map[ n ] ] = idx - vtxoffset;
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

  } // end namespace dgf

} // end namespace Dune
