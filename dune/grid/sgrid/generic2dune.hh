// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH
#define DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH

/** \file
 *  \brief Renumber from the Dune subentity numbering to the SGrid-internal one
 *
 * Once upon a time SGrid used a certain system to number its element subentities.
 * That system worked, but the code that implemented it wasn't documented.
 * So in order to understand that hidden numbering system, a second system was
 * implemented additionally, and it was voted that the new system should be the
 * official Dune numbering system.  For a transitional period code was implemented
 * that computed from one system to the other.  All code was eventually updated
 * to use the new system, except for SGrid.  And that is why the transformation
 * code from Dune numbering to SGrid numbering is still here.
 */

#include <dune/common/static_assert.hh>
#include <dune/common/singleton.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  namespace SGridInternal
  {

    template< unsigned int dim >
    class CubeNumberingTable : public Singleton<CubeNumberingTable<dim> >
    {
    protected:
      using Singleton<CubeNumberingTable<dim> >::instance;

      std::vector<unsigned int> generic2dune_[dim+1];

      static unsigned int dynamicGeneric2dune(unsigned int i, unsigned int codim)
      {
        switch(dim)
        {
        case 3 :
          static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
          return (codim == 2 ? edge[ i ] : i);
        case 4 :
          static unsigned int codim2[ 24 ] =
          { 0, 1, 2, 3, 4, 5, 12, 13, 6, 7, 14, 15,
            8, 9, 16, 17, 20, 21, 10, 11, 18, 19, 22, 23 };
          static unsigned int codim3[ 32 ] =
          { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 24, 25,
            18, 19, 26, 27, 12, 13, 14, 15, 20, 21, 28, 29, 22, 23, 30, 31 };
          if (codim == 2)
            return codim2[i];
          else if (codim == 3)
            return codim3[i];
          else
            return i;
        default :
          return i;
        }
      }

    public:
      CubeNumberingTable ()
      {
        for(unsigned int codim=0; codim<=dim; ++codim)
        {
          generic2dune_[codim].resize(ReferenceElements<double, dim>::cube().size(codim));
          for(unsigned int i=0; i<generic2dune_[codim].size(); ++i)
            generic2dune_[codim][i] = dynamicGeneric2dune(i, codim);
        }
      }

      static unsigned int generic2dune (unsigned int i, unsigned int codim )
      {
        assert(codim <= dim);
        return instance().generic2dune_[codim][i];
      }
    };

  }

}

#endif // #ifndef DUNE_GRID_SGRIDINTERNAL_GENERIC2DUNE_HH
