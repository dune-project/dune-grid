// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/gridparameter.hh>

namespace Dune
{

  namespace dgf
  {

    // GridParameterBlock
    // ------------------

    GridParameterBlock::GridParameterBlock ( std :: istream &in )
      : BasicBlock( in, "GridParameter" ),
        foundFlags_( 0 ),
        name_( "Unnamed Grid" ), // default value (used if name is empty)
        dumpFileName_( "" ),
        markLongestEdge_( false ) // default value
    {
      if( isempty() )
        return;

      // check name
      if( findtoken( "name" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
          name_ = entry;
        else
          dwarn << "GridParameterBlock: Found keyword 'name' without value." << std::endl;
        foundFlags_ |= foundName;
      }

      // get file name
      if ( findtoken( "dumpfilename" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
          dumpFileName_ = entry;
        else
          dwarn << "GridParameterBlock: Found keyword 'dumpFileName' without value." << std::endl;
        foundFlags_ |= foundDumpFileName;
      }

      // check for markLongestEdge
      if( findtoken( "refinementedge" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
        {
          makeupcase( entry );
          if( entry == "LONGEST" )
            markLongestEdge_ = true;
          else if( entry != "ARBITRARY" )
            dwarn << "GridParameterBlock: Invalid value for keyword 'refinementedge': " << entry << std::endl;
        }
        else
          dwarn << "GridParameterBlock: Found keyword 'refinementedge' without value." << std::endl;
        foundFlags_ |= foundLongestEdge;
      }

    }

  } // end namespace dgf

} // end namespace Dune
