// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/general.hh>

namespace Dune
{

  namespace dgf
  {

    // GeneralBlock
    // ---------

    GeneralBlock :: GeneralBlock
      ( std::istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid )
      : BasicBlock( in, "General" ),
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
        for( size_t i = 0, count = 0; i < map.size(); ++i , ++count )
        {
          int x;
          if( !getnextentry( x ) )
          {
            // if we don't find another entry thats it
            map.resize( count );
            break ;
          }
          map[ i ] = x;
        }
      }
    }


    int GeneralBlock :: getDimGrid ()
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


    int GeneralBlock :: get ( std :: vector< std :: vector< unsigned int> > &elements,
                              std :: vector< std :: vector< double > > &params,
                              int &nofp )
    {
      nofp = nofparams;
      reset();

      std :: vector< unsigned int > element( 1 << dimgrid );
      std :: vector< double > param( nofparams );
      int nofcubes = 0;
      for( ; next( element, param ); ++nofcubes )
      {
        elements.push_back( element );
        if( nofparams > 0 )
          params.push_back( param );
      }
      return nofcubes;
    }


    bool GeneralBlock :: next ( std :: vector< unsigned int > &element,
                                std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      for( std :: size_t n = 0; n <element.size(); ++n )
      {
        int idx;
        if( !getnextentry( idx ) )
        {
          // if no entry is avialable and n > 0
          // just resize the element vector
          if( n > 0 )
          {
            element.resize( n ) ;
          }
        }
        if( (vtxoffset > idx) || (idx >= int(nofvtx + vtxoffset)) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid vertex index "
                                  << "(" << idx << " not in ["
                                  << vtxoffset << ", " << (nofvtx + vtxoffset) << "[)" );
        }
        element[ map[ n ] ] = idx - vtxoffset;
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
