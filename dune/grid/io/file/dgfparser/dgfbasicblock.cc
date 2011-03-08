// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfbasicblock.hh>

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


    void BasicBlock :: getblock ( std :: istream &in )
    {
      linecount = 0;
      while( in.good() )
      {
        std :: string line;
        getline( in, line );

        std :: istringstream linestream( line );
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
        std :: string line;
        getline( in, line );

        // strip comments
        if( !line.empty() )
        {
          std::size_t comment = line.find( '%' );
          if( comment != std::string::npos )
            line.erase( comment );
        }
        if( line.empty() )
          continue;

        std :: istringstream linestream( line );
        char test = 0;
        linestream >> test;
        if( test == '#' )
          return;

        ++linecount;
        block << line << std :: endl;
      }
      DUNE_THROW( DGFException, "Error reading from stream." );
    }


    // get next line and store in string stream
    bool BasicBlock :: getnextline ()
    {
      getline( block, oneline );
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
