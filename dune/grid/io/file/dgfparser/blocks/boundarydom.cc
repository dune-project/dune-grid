// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/boundarydom.hh>

namespace Dune
{

  namespace dgf
  {

    // BoundaryDomBlock
    // ----------------

    BoundaryDomBlock :: BoundaryDomBlock(std::istream& in,int cdimworld )
      : BasicBlock( in, "boundarydomain" ),
        dimworld_( cdimworld ),
        counter_( -1 ),
        default_( 0 ),
        ndomains_( 0 )
    {
      if ( !isactive() )
        return;

      assert( cdimworld > 0 );

      if (findtoken("default"))
      {
        int id;
        BoundaryParameter parameter = DGFBoundaryParameter::defaultValue();

        // try to read boundary id
        if ( getnextentry( id ) )
        {
          if( id <= 0 )
          {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      non-positive boundary id (" << id << ") read!");
          }

          // check for parameter
          std::string currentline = line.str();
          std::size_t delimiter = currentline.find( DGFBoundaryParameter::delimiter );
          if( delimiter != std::string::npos )
          {
            parameter =
              DGFBoundaryParameter::convert( currentline.substr( delimiter+1, std::string::npos ) );
          }

          // create default domain data
          default_ = new DomainData( id, parameter, true );

        }
      }

      readBlock();
      reset();
    }

    void BoundaryDomBlock :: readBlock ()
    {
      reset();
      assert( ok() );

      while( getnextline() )
      {
        int id;
        BoundaryParameter parameter = DGFBoundaryParameter::defaultValue();

        if( getnextentry( id ) )
        {
          if( id <= 0 )
          {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      non-positive boundary id (" << id << ") read!");
          }

          // check for parameter
          std::string currentline = line.str();
          std::size_t delimiter = currentline.find( DGFBoundaryParameter::delimiter );
          if( delimiter != std::string::npos )
          {
            parameter =
              DGFBoundaryParameter::convert( currentline.substr( delimiter+1, std::string::npos ) );
          }

          DomainData data( id, parameter );

          // read vertices
          std::vector< double > left( dimworld_ ), right( dimworld_ );
          double x;
          int n = 0;

          while ( getnextentry( x ) )
          {
            if ( 0 <= n && n < dimworld_ )
              left.at( n ) = x;
            else if ( dimworld_ <= n && n < 2*dimworld_ )
            {
              right.at( n-dimworld_ ) = x;
              if ( right.at( n-dimworld_ )< left.at( n-dimworld_ ) )
              {
                DUNE_THROW(DGFException,
                           "ERROR in " << *this
                                       << "      second coordinate smaller than first coordinate: "
                                       << right.at(n-dimworld_)
                                       << " read but expected value larger or equal to "
                                       << left.at(n-dimworld_)
                                       << std::endl << "Line was: '" << line.str() << "'");
              }
            }
            n++;
          }
          bool goodline = ( n == dimworld_*2 );
          if ( !goodline )
          {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      wrong number of coordinates: "
                                   << n << " read but expected 2*" << dimworld_
                                   << std::endl << "Line was: '" << line.str() << "'");
          }

          Domain domain( left, right, data );
          domains_.push_back( domain );

        }
      }

      ndomains_ = domains_.size();
    }

    bool BoundaryDomBlock :: hasParameter () const
    {
      for( int i = 0; i < ndomains_; ++i )
      {
        if( domains_[ i ].data().hasParameter() )
          return true;
      }

      if( hasDefaultData() )
        return defaultData()->hasParameter();

      return false;
    }

  } // end namespace dgf

} // end namespace Dune
