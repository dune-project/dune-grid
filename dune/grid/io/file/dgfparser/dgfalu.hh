// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERALU_HH
#define DUNE_DGFPARSERALU_HH

// only include if ALUGrid is used
#if defined ENABLE_ALUGRID

#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfprojectionblock.hh>

namespace Dune
{

  /** \cond */
  template <int dim,int dimworld>
  class MacroGrid::Impl<ALUCubeGrid<dim,dimworld> > {
    friend class MacroGrid::Impl< ALUConformGrid<dim,dimworld> >;
  public:
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    static ALUCubeGrid<dim,dimworld>*
    generate(MacroGrid& mg,
             const char* filename,
             MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );

    template <class GridImp>
    static GridImp* callDirectly( const char* gridname,
                                  const int rank,
                                  const char *filename,
                                  MPICommunicatorType communicator );

    static bool fileExists ( const char *fileName )
    {
      std :: ifstream testfile( fileName );
      if( !testfile )
        return false;
      testfile.close();
      return true;
    }
  };
  /** \endcond */

  /** \cond */
  template <int dim,int dimworld>
  class MacroGrid::Impl< ALUSimplexGrid<dim,dimworld> >
  {
  public:
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    static ALUSimplexGrid<dim,dimworld>*
    generate(MacroGrid& mg,
             const char* filename,
             MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() )
    {
      // call generic generation of confrom grid (is the same)
      return MacroGrid::Impl< ALUConformGrid<dim,dimworld> > :: template
             generate< ALUSimplexGrid<dim,dimworld> >
               ("ALUSimplexGrid< 2 , 2 >", mg, filename, MPICOMM );
    }
  };
  /** \endcond */


  /** \cond */
  template <int dim,int dimworld>
  class MacroGrid::Impl< ALUConformGrid<dim,dimworld> >
  {
    friend class MacroGrid::Impl< ALUSimplexGrid<dim,dimworld> >;
  public:
    // type of MPI communicator
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    static ALUConformGrid<dim,dimworld>*
    generate(MacroGrid& mg,
             const char* filename,
             MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() )
    {
      return MacroGrid::Impl< ALUConformGrid<dim,dimworld> > ::
             template generate < ALUConformGrid<dim,dimworld> >
               ("ALUConformGrid< 2 , 2 >", mg, filename, MPICOMM );
    }

    // generic 2d grid generation method
    template <class GridImp>
    static GridImp*
    generate(const char* gridname,
             MacroGrid& mg,
             const char* filename,
             MPICommunicatorType MPICOMM );


    static bool fileExists ( const char *fileName )
    {
      return MacroGrid::Impl< ALUCubeGrid<dim,dimworld> > :: fileExists( fileName );
    }
  };
  /** \endcond */

  // DGFGridInfo (specialization for ALUGrid)
  // ----------------------------------------

  /** \cond */
  template<>
  struct DGFGridInfo< ALUCubeGrid< 3, 3 > >
  {
    static int refineStepsForHalf () { return 1; }
    static double refineWeight () { return 0.125; }
  };

  template<>
  struct DGFGridInfo< ALUSimplexGrid< 3, 3 > >
  {
    static int refineStepsForHalf () { return 1; }
    static double refineWeight () { return 0.125; }
  };

  template<>
  struct DGFGridInfo< ALUSimplexGrid< 2, 2 > >
  {
    static int refineStepsForHalf () { return 1; }
    static double refineWeight () { return 0.25; }
  };

  template<>
  struct DGFGridInfo< Dune::ALUConformGrid< 2, 2 > >
  {
    static int refineStepsForHalf () { return 2; }
    static double refineWeight () { return 0.5; }
  };
  /** \endcond */

}

#include "dgfalu.cc"

#endif // #if defined ENABLE_ALUGRID

#endif // #ifndef DUNE_DGFPARSERALU_HH
