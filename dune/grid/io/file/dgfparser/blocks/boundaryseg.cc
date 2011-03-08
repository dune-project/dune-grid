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

    BoundarySegBlock :: BoundarySegBlock ( std :: istream &in, int pnofvtx,
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
    :: get( std :: map< DGFEntityKey< unsigned int>, int > &facemap,
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

          DGFEntityKey< unsigned int > key( bound, false );

          facemap[key] = bndid;
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

                DGFEntityKey< unsigned int > key( bound, false );
                facemap[key] = bndid;
                ++lnofbound;
              }
            }
            else
            {
              for (int i=2; i<=size(); i++)
              {
                bound[0] = p[i-1];
                bound[1] = p[i%size()];
                DGFEntityKey< unsigned int > key( bound, false );
                facemap[key] = bndid;
                ++lnofbound;
              }
            }
          }
        }
        else {
          if (dimworld==3) {
            DGFEntityKey< unsigned int > key( p, false );
            facemap[key] = bndid;
            ++lnofbound;
          } else {
            std :: vector< unsigned int > k( 2 );
            for (size_t i=0; i<p.size()-1; i++) {
              k[0]=p[i];
              k[1]=p[(i+1)%p.size()];
              DGFEntityKey< unsigned int > key( k, false );
              facemap[key] = bndid;
              ++lnofbound;
            }
          }
        }
      }
      return lnofbound;
    }


    bool BoundarySegBlock :: next ()
    {
      assert(ok());
      int n=0;
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      p.clear();
      int x;
      if (getnextentry(x)) {
        bndid = x;
        if (bndid<=0)
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      non-positive boundary id (" << bndid << ") read!");
        }
        while (getnextentry(x))
        {
          p.push_back(x);
          n++;
        }
        // goodline=(n==dimworld+1);
        goodline=true;
        return goodline;
      } else return next();
    }

  } // end namespace dgf

} // end namespace Dune
