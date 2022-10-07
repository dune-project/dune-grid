// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {

    // BasicBlock
    // ----------

    BasicBlock :: BasicBlock ( std :: istream &in, const char *id )
      : pos(-1),
        active(false),
        empty(true),
        identifier(id),
        linecount(0)
    {
      makeupcase( identifier );
      in.clear();
      in.seekg(0);
      if (!in) {
        DUNE_THROW(DGFException,
                   "file not found in BasicBlock::BasicBlock");
      }
      getblock( in );
      empty = (linecount == 0);
      if (active && !empty) {
        //linecount=countlines();
        reset();
      }
      in.clear();
      in.seekg(0);
    }

    // read the current block which is ended by a line starting
    // with a # symbol.
    void BasicBlock :: getblock ( std :: istream &in )
    {
      linecount = 0;
      while( in.good() )
      {
        std::string curLine;
        getline( in, curLine );

        std::istringstream linestream( curLine );
        std :: string id;
        linestream >> id;

        makeupcase( id );
        if( id == identifier )
          break;
      }
      if( in.eof() )
        return;

      active = true;
      while( in.good() )
      {
        std::string curLine;
        getline( in, curLine );

        // strip comments
        if( !curLine.empty() )
        {
          std::size_t comment = curLine.find( '%' );
          if( comment != std::string::npos )
            curLine.erase( comment );
        }
        if( curLine.empty() )
          continue;

        std :: istringstream linestream( curLine );
        char test = 0;
        linestream >> test;
        if( test == '#' )
          return;

        ++linecount;
        block_ << curLine << std :: endl;
      }
      DUNE_THROW( DGFException,
                  "Error reading from stream, expected \"#\" to end the block." );
    }


    // get next line and store in string stream
    bool BasicBlock :: getnextline ()
    {
      getline( block_, oneline );
      line.clear();
      line.str( oneline );
      ++pos;
      return !oneline.empty();
    }


    bool BasicBlock :: gettokenparam ( std :: string token, std :: string &entry )
    {
      reset();
      makeupcase( token );
      while( getnextline() )
      {
        std :: string ltoken;
        line >> ltoken;
        makeupcase( ltoken );
        if( ltoken == token )
        {
          getline( line, entry );
          return true;
        }
      }
      return false;
    }


    bool BasicBlock :: findtoken ( std :: string token )
    {
      reset();
      makeupcase( token );
      while( getnextline() )
      {
        std :: string ltoken;
        line >> ltoken;
        makeupcase( ltoken );
        if( ltoken == token )
          return true;
      }
      return false;
    }

  } // end namespace dgf

} // end namespace Dune
