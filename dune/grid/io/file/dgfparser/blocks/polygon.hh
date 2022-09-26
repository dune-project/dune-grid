// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYGON_HH
#define DUNE_POLYGON_HH

#include <iostream>
#include <vector>

#include <dune/common/typetraits.hh>
#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {

    // PolygonBlock
    // ------------

    struct PolygonBlock
      : public BasicBlock
    {
      PolygonBlock ( std::istream &in, int numVtx, int vtxOfs )
        : BasicBlock( in, "Polygon" ), vtxBegin_( vtxOfs ), vtxEnd_( vtxOfs + numVtx )
      {}

      int get ( std::vector< std::vector< int > > &polygons )
      {
        reset();
        std::vector< int > polygon;
        while( getnextline() )
        {
          polygon.clear();
          for( int vtxIdx; getnextentry( vtxIdx ); )
          {
            if( (vtxBegin_ > vtxIdx) || (vtxIdx >= vtxEnd_) )
              DUNE_THROW( DGFException, "Error in " << *this << ": Invalid vertex index (" << vtxIdx << " not int [" << vtxBegin_ << ", " << vtxEnd_ << "[)" );
            polygon.push_back( vtxIdx - vtxBegin_ );
          }

          polygons.push_back( polygon );
        }
        return polygons.size();
      }

    protected:
      int vtxBegin_, vtxEnd_;
    };

  } // namespace dgf
} // end namespace Dune

#endif // #ifndef DUNE_POLYGON_HH
