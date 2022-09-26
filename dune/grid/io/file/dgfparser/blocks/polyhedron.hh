// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRON_HH
#define DUNE_POLYHEDRON_HH

#include <algorithm>

#include <dune/grid/io/file/dgfparser/blocks/polygon.hh>

namespace Dune
{

  namespace dgf
  {

    // PolyhedronBlock
    // ---------------

    struct PolyhedronBlock
      : public BasicBlock
    {
      explicit PolyhedronBlock ( std::istream &in, int numPolys )
        : BasicBlock( in, "Polyhedron" ), numPolys_( numPolys )
      {}

      int get ( std::vector< std::vector< int > > &polyhedra )
      {
        reset();
        std::vector< int > polyhedron;
        int minPolyId = 1;
        while( getnextline() )
        {
          polyhedron.clear();
          for( int polyIdx; getnextentry( polyIdx ); )
          {
            if( (polyIdx < 0) || (polyIdx > numPolys_) )
              DUNE_THROW( DGFException, "Error in " << *this << ": Invalid polygon index (" << polyIdx << " not int [0, " << numPolys_ << "])" );

            minPolyId = std::min( minPolyId, polyIdx );
            polyhedron.push_back( polyIdx );
          }

          polyhedra.push_back( polyhedron );
        }

        // subtract minimal number to have 0 starting numbering
        if( minPolyId > 0 )
        {
          const size_t polySize = polyhedra.size();
          for( size_t i=0; i<polySize; ++i )
          {
            const size_t pSize = polyhedra[ i ].size();
            for( size_t j=0; j<pSize; ++j )
            {
              polyhedra[ i ][ j ] -= minPolyId;
            }
          }
        }
        return polyhedra.size();
      }

    protected:
      const int numPolys_;
    };

  } // namespace dgf
} // end namespace Dune

#endif // #ifndef DUNE_POLYHEDRON_HH
