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

  // DGFGridFactory for AluSimplexGrid
  // ------------------------------

  // template< int dim, int dimworld > // for a first version
  template < int dim >
  struct DGFGridFactory< ALUSimplexGrid< dim, 3 > >
  {
    typedef ALUSimplexGrid<3,3>  Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;
    typedef Dune::GridFactory<Grid> GridFactory;

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : factory_(comm),
        dgf_( 0, 1 )
    {
      generate( filename, comm );
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return factory_.wasInserted( intersection );
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      if( codim == 0 )
        return dgf_.nofelparams;
      else if( codim == dimension )
        return dgf_.nofvtxparams;
      else
        return 0;
    }

    std::vector< double > &parameter ( const Element &element )
    {
      if( numParameters< 0 >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.elParams[ factory_.insertionIndex( element ) ];
    }

    std::vector< double > &parameter ( const Vertex &vertex )
    {
      if( numParameters< dimension >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.vtxParams[ factory_.insertionIndex( vertex ) ];
    }

  private:
    void generate( const std::string &filename,
                   MPICommunicatorType comm );

    Grid *grid_;
    GridFactory factory_;
    DuneGridFormatParser dgf_;
  };
}

#include "dgfalu.cc"

#endif // #if defined ENABLE_ALUGRID

#endif // #ifndef DUNE_DGFPARSERALU_HH
