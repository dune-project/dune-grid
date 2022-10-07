// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_BASICBLOCK_HH
#define DUNE_DGF_BASICBLOCK_HH

#include <cassert>
#include <cctype>
#include <iostream>
#include <string>
#include <sstream>

#include <dune/common/stdstreams.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

namespace Dune
{

  namespace dgf
  {

    inline void makeupcase( std :: string &s )
    {
      for (size_t i=0; i<s.size(); i++)
        s[i]=std::toupper(s[i]);
    }

    class BasicBlock
    {
      int pos;                   // line number
      bool active;               // block was found
      bool empty;                // block was found but was empty
      std::string identifier;    // identifier of this block
      int linecount;             // total number of lines in the block
      std::stringstream block_;  // the block itself
      std::string oneline;       // the active line in the block

      // get the block (if it exists)
      void getblock ( std::istream &in );

      // count the number of lines in the block
      // int countlines ();

    protected:
      std::stringstream line; // the active line as string buffer
                              // for use in the derived classes

      // go back to beginning of block
      void reset ()
      {
        pos = -1;
        block_.clear();
        block_.seekg( 0 );
      }

      // get next line and store in string stream
      bool getnextline ();

      // get next entry in line
      template< class ENTRY >
      bool getnextentry( ENTRY &entry )
      {
        line >> entry;
        return static_cast< bool >( line );
      }

      bool gettokenparam ( std :: string token, std :: string &entry );
      bool findtoken( std :: string token );

    public:
      // search for block in file and store in buffer
      BasicBlock ( std::istream &in, const char* id );

      // some information on this block
      bool isactive ()
      {
        return active;
      }

      bool isempty ()
      {
        return empty;
      }

      int &noflines ()
      {
        return linecount;
      }

      int linenumber ()
      {
        return pos;
      }

      const std::string & id () const
      {
        return identifier;
      }

      // for error messages
      friend std :: ostream &operator<< ( std :: ostream &os, const BasicBlock &b )
      {
        return os << "block " << b.identifier << " (line " << b.pos << ")";
      }

    };

  } // end namespace dgf

} // end namespace Dune

#endif
