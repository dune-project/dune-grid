// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/interval.hh>

namespace Dune
{

  namespace dgf
  {

    // IntervalBlock
    // -------------

    IntervalBlock::IntervalBlock ( std::istream &in )
      : BasicBlock( in, "Interval" ),
        intervals_( 0 ),
        map_(),
        good_( false ),
        dimw_( 0 )
    {
      if( !isactive() )
        return;

      const bool hasMapping = findtoken("map");

      if( hasMapping )
      {
        int x;
        map_.clear();
        while( getnextentry( x ) )
        {
          map_.push_back( x );
        }
      }
      else
      {
        reset();
      }

      getnextline();

      for( double x; getnextentry( x ); ++dimw_ ) ;
      if( dimw_ == 0 )
      {
        DUNE_THROW( DGFException,
                    "Too few coordinates for point p0 in IntervalBlock" );
      }

      // reset stream
      reset();
      // skip first line if mapping was found
      if( hasMapping )
        getnextline();

      while( next() ) ;

      // set default map
      if( map_.empty() )
      {
        map_.resize( dimw_ );
        for( int i=0; i<dimw_; ++i )
          map_[ i ] = i;
      }

      // check mapping
      if( int(map_.size()) != dimw_ )
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": "
                                << "Incomplete ijk mapping "
                                << "(got " << map_.size() << " entries, "
                                << "expected " << dimw_ << " entries." );
      for( int i=0; i<dimw_; ++i )
      {
        if( map_[ i ] < 0 || map_[ i ] >= dimw_ )
          DUNE_THROW( DGFException, "Error in " << *this << ": ijk mapping is not a permutation of {0,...," << dimw_-1 << "}");
      }
    }


    int IntervalBlock::getVtx ( int block, std::vector< std::vector< double > > &vtx ) const
    {
      dverb << "reading vertices for interval " << block << "... ";

      Interval interval( get( block ), map_ );

      size_t old_size = vtx.size();
      vtx.resize( old_size + nofvtx( block ) );
      for( size_t i = old_size; i < vtx.size(); ++i )
        vtx[ i ].resize( dimw() );

      size_t m = old_size;
      std::vector< int > i( dimw() );
      const int end = dimw()-1;
      int k = end;
      for( i[ end ] = 0; i[ end ] <= interval.n[ end ]; )
      {
        // go all the way down
        for( ; k > 0; --k )
          i[ k-1 ] = 0;

        assert( m < vtx.size() );
        for( int j = 0; j < dimw(); ++j ) {
          // apply mapping to get correct coordinates
          vtx[ m ][ map_[ j ] ] = interval.p[ 0 ][ j ] + double(i[ j ])*interval.h[ j ];
        }
        ++m;

        // increase i[ k ] and go up for all finished loops
        for( ; (++i[ k ] > interval.n[ k ]) && (k < end); ++k ) ;
      }
      assert( m == vtx.size() );

      dverb << "[done]" << std::endl;
      return m - old_size;;
    }


    int IntervalBlock::getHexa ( int block, std::vector< std::vector< unsigned int > > &cubes, int offset ) const
    {
      dverb << "generating cubes for interval " << block << "... ";

      Interval interval( get( block ), map_ );
      const int verticesPerCube = 1 << dimw();

      size_t old_size = cubes.size();
      cubes.resize( old_size + nofhexa( block ) );
      for( size_t i = old_size; i < cubes.size(); ++i )
        cubes[ i ].resize( verticesPerCube );

      size_t m = old_size;
      std::vector< int > i( dimw() );
      const int end = dimw()-1;
      int k = end;
      for( i[ end ] = 0; i[ end ] < interval.n[ end ]; )
      {
        // go all the way down
        for( ; k > 0; --k )
          i[ k-1 ] = 0;

        assert( m < cubes.size() );
        for( int j = 0; j < verticesPerCube; ++j )
        {
          cubes[ m ][ j ] = offset;
          int factor = 1;
          for( int d = 0; d < dimw(); ++d )
          {
            cubes[ m ][ j ] += factor*(i[ d ] + ((j >> d) & 1));
            factor *= interval.n[ d ]+1;
          }
        }
        ++m;

        // increase i[ k ] and go up for all finished loops
        for( ; (++i[ k ] >= interval.n[ k ]) && (k < end); ++k ) ;
      }
      assert( m == cubes.size() );

      dverb << "[done]" << std::endl;
      return m - old_size;
    }


    template< class T >
    void IntervalBlock::parseLine ( std::vector< T > &v )
    {
      getnextline();
      v.resize( dimw_ );
      for( int i = 0; i < dimw_; ++i )
      {
        if( !getnextentry( v[ i ] ) )
          DUNE_THROW( DGFException, "ERROR in " << *this << ": Not enough values." );
      }
    }


    bool IntervalBlock::next ()
    {
      if (linenumber()==noflines()-1)
      {
        good_=false;
        return good_;
      }

      Interval interval;
      parseLine( interval.p[ 0 ] );
      parseLine( interval.p[ 1 ] );
      parseLine( interval.n );

      //find real upper and lower edge and calculate cell width
      interval.h.resize( dimw_ );
      for( int i = 0; i < dimw_; ++i )
      {
        double &left = interval.p[ 0 ][ i ];
        double &right = interval.p[ 1 ][ i ];
        const int &n = interval.n[ i ];

        if( left > right )
        {
          double dummy = left;
          left = right;
          right = dummy;
        }

        interval.h[ i ] = (right - left) / double( n );
        assert( interval.h[ i ] > 0);
      }
      intervals_.push_back( interval );

      dverb << interval << std::endl;

      good_ = true;
      return good_;
    }

  } // end namespace dgf

} // end namespace Dune
