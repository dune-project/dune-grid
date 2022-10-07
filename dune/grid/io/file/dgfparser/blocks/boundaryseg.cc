// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/boundaryseg.hh>

namespace Dune
{

  namespace dgf
  {

    // BoundarySegBlock
    // ----------------

    BoundarySegBlock :: BoundarySegBlock ( std :: istream &in, int /* pnofvtx */,
                                           int pdimworld,bool psimplexgrid )
      : BasicBlock(in,"boundarysegments"),
        dimworld(pdimworld),
        goodline(true),
        p(),
        bndid(-1),
        simplexgrid(psimplexgrid)
    {
      if (!isactive())
        return;
      assert(dimworld>0);
      next();
    }


    int BoundarySegBlock
    :: get( std :: map< EntityKey, BndParam > &facemap,
            bool fixedsize,
            int vtxoffset )
    {
      static int cube2simplex[3][3] = {
        {0,1,3},
        {0,2,3},
        {1,2,3}
      };

      int lnofbound;
      int face = ElementFaceUtil::faceSize(dimworld,simplexgrid);
      for (lnofbound=0; ok(); next())
      {
        for (size_t i=0; i<p.size(); i++)
        {
          p[i] -= vtxoffset;
        }
        if (fixedsize)
        {
          if ((dimworld==2 && size()< 2) ||
              (dimworld==3 && simplexgrid && size()!=3 && size()!=4) ||
              (dimworld==3 && !simplexgrid && size()!=4))
            continue;
          std :: vector< unsigned int > bound( face );
          for (int j=0; j<face; j++) {
            bound[j] = p[j];
          }

          EntityKey key( bound, false );

          facemap[key] = BndParam( bndid, parameter );
          ++lnofbound;

          if (size()>face)
          {
            assert(dimworld==2 || face==3);
            if (dimworld==3) {
              for (int i=0; i<3; i++)
              {
                for (int j=0; j<face; j++)
                {
                  bound[j] = p[cube2simplex[i][j]];
                }

                EntityKey helperKey( bound, false );
                facemap[helperKey] = BndParam( bndid, parameter );
                ++lnofbound;
              }
            }
            else
            {
              for (int i=2; i<=size(); i++)
              {
                bound[0] = p[i-1];
                bound[1] = p[i%size()];
                EntityKey helperKey( bound, false );
                facemap[helperKey] = BndParam( bndid, parameter );
                ++lnofbound;
              }
            }
          }
        }
        else {
          if (dimworld==3) {
            EntityKey key( p, false );
            facemap[key] = BndParam( bndid, parameter );
            ++lnofbound;
          } else {
            std :: vector< unsigned int > k( 2 );
            for (size_t i=0; i<p.size()-1; i++) {
              k[0]=p[i];
              k[1]=p[(i+1)%p.size()];
              EntityKey key( k, false );
              facemap[key] = BndParam( bndid, parameter );
              ++lnofbound;
            }
          }
        }
      }
      return lnofbound;
    }


    bool BoundarySegBlock :: next ()
    {
      assert( ok() );
      getnextline();

      if (linenumber()==noflines())
      {
        goodline=false;
        return goodline;
      }

      // clear data
      p.clear();
      parameter = DGFBoundaryParameter::defaultValue();

      // get active line
      std::string currentline = line.str();
      if( !currentline.empty() )
      {
        // find delimiter and split line
        std::size_t delimiter = currentline.find( DGFBoundaryParameter::delimiter );
        std::string left = currentline.substr( 0, delimiter );

        // read boundary id and boundary vertices from left hand side
        {
          std::istringstream lstream( left );
          assert( !left.empty() );

          int x;
          lstream >> x;
          bndid = x;
          if ( bndid <= 0 )
          {
            DUNE_THROW( DGFException,
                        "ERROR in " << *this
                                    << "      non-positive boundary id (" << bndid << ") read!" );
          }

          while( lstream >> x )
            p.push_back( x );
        }

        // read parameter from right hand side
        if( delimiter != std::string::npos )
        {
          std::string strParam = currentline.substr( delimiter+1, std::string::npos );
          parameter = DGFBoundaryParameter::convert( strParam );
        }

        goodline=true;
        return goodline;
      }
      else
        return next();
    }

  } // end namespace dgf

} // end namespace Dune
