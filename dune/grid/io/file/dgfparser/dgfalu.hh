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
  template < class G >
  struct DGFBaseFactory
  {
    typedef G Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;
    typedef Dune::GridFactory<Grid> GridFactory;

    DGFBaseFactory ()
      : factory_( ),
        dgf_( 0, 1 )
    {}

    explicit DGFBaseFactory ( MPICommunicatorType comm )
      : factory_(comm),
        dgf_( rank(comm), size(comm) )
    {}

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

  protected:
    static Grid* callDirectly( const char* gridname,
                               const int rank,
                               const char *filename,
                               MPICommunicatorType communicator )
    {
  #if ALU3DGRID_PARALLEL
      std :: stringstream tmps;
      tmps << filename << "." << rank;
      const std :: string &tmp = tmps.str();

      if( fileExists( tmp.c_str() ) )
        return new Grid( tmp.c_str(), communicator );
  #endif
      if( fileExists( filename ) )
      {
        if( rank == 0 )
          return new Grid( filename );
        else
          return new Grid( );
      }
      DUNE_THROW( GridError, "Unable to create " << gridname << " from '"
                                                 << filename << "'." );
    }
    static bool fileExists ( const char *fileName )
    {
      std :: ifstream testfile( fileName );
      if( !testfile )
        return false;
      testfile.close();
      return true;
    }
    static int rank( MPICommunicatorType MPICOMM )
    {
      int rank = 0;
#if HAVE_MPI
      MPI_Comm_rank( MPICOMM, &rank );
#endif
      return rank;
    }
    static int size( MPICommunicatorType MPICOMM )
    {
      int size = 1;
#if HAVE_MPI
      MPI_Comm_size( MPICOMM, &size );
#endif
      return size;
    }
    Grid *grid_;
    GridFactory factory_;
    DuneGridFormatParser dgf_;
  };

  template <>
  struct DGFGridFactory< ALUSimplexGrid<3,3> > :
    public DGFBaseFactory< ALUSimplexGrid<3,3> >
  {
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< ALUSimplexGrid<3,3> >( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUSimplexGrid<3,3> >( comm )
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
        grid_ = callDirectly( "ALUSimplexGrid< 3 , 3 >", rank( comm ), filename.c_str(), comm );
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  template <>
  struct DGFGridFactory< ALUCubeGrid<3,3> > :
    public DGFBaseFactory< ALUCubeGrid<3,3> >
  {
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< ALUCubeGrid<3,3> >( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUCubeGrid<3,3> >( comm )
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
        grid_ = callDirectly( "ALUCubeGrid< 3 , 3 >", rank( comm ), filename.c_str(), comm );
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  template <>
  struct DGFGridFactory< ALUConformGrid<2,2> > :
    public DGFBaseFactory< ALUConformGrid<2,2> >
  {
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< ALUConformGrid<2,2> >()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUConformGrid<2,2> >()
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( fileExists( filename.c_str() ) )
          grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }
  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  template <>
  struct DGFGridFactory< ALUSimplexGrid<2,2> > :
    public DGFBaseFactory< ALUSimplexGrid<2,2> >
  {
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< ALUSimplexGrid<2,2> >()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW(DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUSimplexGrid<2,2> >()
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( fileExists( filename.c_str() ) )
          grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

}

#include "dgfalu.cc"

#endif // #if defined ENABLE_ALUGRID

#endif // #ifndef DUNE_DGFPARSERALU_HH
