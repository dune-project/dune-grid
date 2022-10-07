// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_GRIDPARAMETERBLOCK_HH
#define DUNE_DGF_GRIDPARAMETERBLOCK_HH

#include <iostream>
#include <string>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {
    /** \brief Common Grid parameters
        \ingroup DGFGridParameter

         For each grid implementation there is a set of parameters
         that can be passed via the GridParameter block to the momment of
         grid construction.
         Currently implemented common parameters are: \n\n
         1. \b name: The name of the grid ( later returned by the method grid.name() ). \n
         2. \b refinementedge: parameter to specify the refinement edge in simplices.
            Valid values are \b arbitrary (which is the default value) and \b longest
            which marks the longest edge/face of each simplex to be the refinement edge. \n
        See also the \b examplegrid5.dgf file for examples.
     */

    class GridParameterBlock
      : public BasicBlock
    {
    public:
      typedef unsigned int Flags;

      static const Flags foundName = 1 << 0;
      static const Flags foundDumpFileName = 1 << 1;
      static const Flags foundLongestEdge = 1 << 5;

    protected:
      Flags foundFlags_; // supportFlags, this block was created with
      std::string name_; // name of the grid
      std::string dumpFileName_; // name of the grid
      bool markLongestEdge_; // Mark longest edge for AlbertaGrid

    private:
      // copy not implemented
      GridParameterBlock(const GridParameterBlock&);

    public:
      //! constructor: read commmon parameters
      GridParameterBlock ( std::istream &in );

      //! return the name of the grid
      const std::string &name ( const std::string &defaultValue ) const
      {
        if( (foundFlags_ & foundName) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'name' not specified, "
                << "defaulting to '" << defaultValue << "'." << std::endl;
          return defaultValue;
        }
        else
          return name_;
      }

      const std::string &dumpFileName ( ) const
      {
        if( (foundFlags_ & foundDumpFileName) != 0 )
        {
          dwarn << "GridParameterBlock: found Parameter 'dumpfilename', "
                << "dumping file to `" << dumpFileName_ << "'" << std::endl;
        }
        return dumpFileName_;
      }

      //! returns true if longest edge should be marked for AlbertaGrid
      bool markLongestEdge () const
      {
        if( (foundFlags_ & foundLongestEdge) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'refinementedge' not specified, "
                << "defaulting to 'ARBITRARY'." << std::endl;
        }
        return markLongestEdge_;
      }

      // some information
      bool ok()
      {
        return true;
      }
    };


  } // end namespace dgf

} // end namespace Dune

#endif
