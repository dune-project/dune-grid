// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_INTERVALBLOCK_HH
#define DUNE_DGF_INTERVALBLOCK_HH

#include <iostream>
#include <vector>
#include <array>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {

    struct IntervalBlock
      : public BasicBlock
    {
      struct Interval
      {
        Interval() {}
        Interval( const Interval& interval, const std::vector<int>& map )
        {
          copy( interval, map );
        }
        void copy(const Interval& interval, const std::vector<int>& map )
        {
          const int size = map.size();
          p[0].resize( size );
          p[1].resize( size );
          n.resize( size );
          h.resize( size );
          assert( size == int(interval.n.size()) );
          for( int i=0; i<size; ++i )
          {
            p[ 0 ][ i ] = interval.p[ 0 ][ map[ i ] ];
            p[ 1 ][ i ] = interval.p[ 1 ][ map[ i ] ];
            n[ i ] = interval.n[ map[ i ] ];
            h[ i ] = interval.h[ map[ i ] ];
          }
        }
        std::array< std::vector< double >, 2 > p; // lower and upper boundary points
        std::vector< double > h;             // width of the cells in each direction
        std::vector< int > n;                // number of cells in each direction
      };

    private:
      std::vector< Interval > intervals_;
      std::vector< int > map_;
      bool good_;                      //data read correctly
      int dimw_;                       //dimension of world

    public:
      explicit IntervalBlock ( std::istream &in );

      void get ( std::vector< std::vector< double > > &vtx, int &nofvtx,
                 std::vector< std::vector< unsigned int > > &simplex, int &nofsimpl )
      {
        for( size_t i = 0; i < intervals_.size(); ++i )
        {
          int oldvtx = nofvtx;
          nofvtx += getVtx( i, vtx );
          nofsimpl += getHexa( i, simplex, oldvtx );
        }
      }

      void get ( std::vector< std::vector< double > > &vtx, int &nofvtx )
      {
        for( size_t i = 0; i < intervals_.size(); ++i )
          nofvtx += getVtx( i, vtx );
      }

      const Interval &get ( int block ) const
      {
        return intervals_[ block ];
      }

      int numIntervals () const
      {
        return intervals_.size();
      }

      int dimw () const
      {
        return dimw_;
      }

      int getVtx ( int block, std::vector< std::vector< double > > &vtx ) const;
      int getHexa ( int block, std::vector< std::vector< unsigned int > > &cubes,
                    int offset = 0 ) const;

      int nofvtx ( int block ) const
      {
        const Interval &interval = get( block );
        int n = 1;
        for( int i = 0; i < dimw_; ++i )
          n *= (interval.n[ i ] + 1);
        return n;
      }

      int nofhexa ( int block ) const
      {
        const Interval &interval = get( block );
        int n = 1;
        for( int i = 0; i < dimw_; ++i )
          n *= interval.n[ i ];
        return n;
      }

    private:
      template< class T >
      void parseLine ( std::vector< T > &v );

      bool next ();
    };

    inline std::ostream &
    operator<< ( std::ostream &out, const IntervalBlock::Interval &interval )
    {
      if( interval.p[ 0 ].empty() || interval.p[ 1 ].empty() || interval.n.empty() )
        return out << "Interval {}";

      out << "Interval { p0 = (" << interval.p[ 0 ][ 0 ];
      for( size_t i = 1; i < interval.p[ 0 ].size(); ++i )
        out << ", " << interval.p[ 0 ][ i ];
      out << "), p1 = (" << interval.p[ 1 ][ 0 ];
      for( size_t i = 1; i < interval.p[ 1 ].size(); ++i )
        out << ", " << interval.p[ 1 ][ i ];
      out << "), n = (" << interval.n[ 0 ];
      for( size_t i = 1; i < interval.n.size(); ++i )
        out << ", " << interval.n[ i ];
      return out << ") }";
    }

  } // end namespace dgf

} // end namespace Dune

#endif
