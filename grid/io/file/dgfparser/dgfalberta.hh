// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERALBERTA_HH
#define DUNE_DGFPARSERALBERTA_HH

#if defined HAVE_ALBERTA

#include <dune/grid/albertagrid.hh>
#include "dgfparser.hh"

namespace Dune
{

  template< int dim, int dimworld >
  class MacroGrid :: Impl< AlbertaGrid< dim, dimworld > >
  {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  public:
    static AlbertaGrid< dim, dimworld > *
    generate( MacroGrid &mg,
              const char* filename,
              MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };


  template< int dim, int dimworld >
  inline AlbertaGrid< dim, dimworld > *
  MacroGrid :: Impl< AlbertaGrid< dim, dimworld > >
  :: generate( MacroGrid &mg, const char* filename, MPICommunicatorType )
  {
    mg.element=Simplex;
    std :: ifstream gridin( filename );

    std :: string str( filename );

    if( mg.readDuneGrid( gridin, dim, dimworld ) )
    {
      /*
         if( mg.dimw != dimworld)
         {
         DUNE_THROW(DGFException,
                   "Macrofile " << filename << " is for dimension " << mg.dimw
                   << " and connot be used to initialize an "
                   << "AlbertaGrid of dimension "
                   << dimworld);
         }
       */
      mg.setOrientation( 0, 1 );
      mg.setRefinement( 0, 1, -1, -1 );
      str += ".albertagrid";
      std :: ofstream out( str.c_str() );
      mg.writeAlberta( out );
    }
    return new AlbertaGrid< dim, dimworld >( str.c_str() );
  }



  template< int dim, int dimworld >
  struct DGFGridInfo< AlbertaGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
      return dim;
    }

    static double refineWeight ()
    {
      return 0.5;
    }
  };

}

#endif

#endif
