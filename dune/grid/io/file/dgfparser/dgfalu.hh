// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERALU_HH
#define DUNE_DGFPARSERALU_HH

// only include if ALUGrid is used
#if defined ENABLE_ALUGRID

#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>
#include <dune/grid/io/file/dgfparser/blocks/projection.hh>

#include <dune/grid/common/intersection.hh>

namespace Dune
{

  // forward declaration
  // -------------------

  template< class GridImp, template< class > class IntersectionImp >
  class Intersection;

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

  template<int dimw>
  struct DGFGridInfo< ALUSimplexGrid< 2, dimw > >
  {
    static int refineStepsForHalf () { return 1; }
    static double refineWeight () { return 0.25; }
  };

  template<int dimw>
  struct DGFGridInfo< ALUCubeGrid< 2, dimw > >
  {
    static int refineStepsForHalf () { return 1; }
    static double refineWeight () { return 0.25; }
  };

  template<int dimw>
  struct DGFGridInfo< Dune::ALUConformGrid< 2, dimw > >
  {
    static int refineStepsForHalf () { return 2; }
    static double refineWeight () { return 0.5; }
  };

  template<int dimg, int dimw, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
  struct DGFGridInfo< Dune::ALUGrid< dimg, dimw, eltype, refinementtype, Comm > >
  {
    static int refineStepsForHalf () { return ( refinementtype == conforming ) ? dimg : 1; }
    static double refineWeight () { return ( refinementtype == conforming ) ? 0.5 : 1.0/(std::pow( 2.0, double(dimg))); }
  };
  /** \endcond */

  // DGFGridFactory for AluGrid
  // --------------------------

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
      : factory_(),
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

    template < class GG, template < class > class II >
    int boundaryId ( const Intersection< GG, II > & intersection ) const
    {
      typedef Dune::Intersection< GG, II > Intersection;
      typename Intersection::EntityPointer inside = intersection.inside();
      const typename Intersection::Entity & entity = *inside;
      const int face = intersection.indexInInside();

      const GenericReferenceElement< void, dimension > &refElement
        = GenericReferenceElements< void, dimension >::general( entity.type() );
      int corners = refElement.size( face, 1, dimension );
      std :: vector< unsigned int > bound( corners );
      for( int i=0; i < corners; ++i )
      {
        const int k =  refElement.subEntity( face, 1, i, dimension );
        bound[ i ] = factory_.insertionIndex( *entity.template subEntity< dimension >( k ) );
      }

      DuneGridFormatParser::facemap_t::key_type key( bound, false );
      const DuneGridFormatParser::facemap_t::const_iterator pos = dgf_.facemap.find( key );
      if( pos != dgf_.facemap.end() )
        return dgf_.facemap.find( key )->second.first;
      else
        return (intersection.boundary() ? 1 : 0);
    }

    template < class GG, template < class > class II >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection< GG, II > & intersection ) const
    {
      typedef Dune::Intersection< GG, II > Intersection;
      typename Intersection::EntityPointer inside = intersection.inside();
      const typename Intersection::Entity & entity = *inside;
      const int face = intersection.indexInInside();

      const GenericReferenceElement< void, dimension > & refElement =
        GenericReferenceElements< void, dimension >::general( entity.type() );
      int corners = refElement.size( face, 1, dimension );
      std :: vector< unsigned int > bound( corners );
      for( int i=0; i < corners; ++i )
      {
        const int k =  refElement.subEntity( face, 1, i, dimension );
        bound[ i ] = factory_.insertionIndex( *entity.template subEntity< dimension >( k ) );
      }

      DuneGridFormatParser::facemap_t::key_type key( bound, false );
      const DuneGridFormatParser::facemap_t::const_iterator pos = dgf_.facemap.find( key );
      if( pos != dgf_.facemap.end() )
        return dgf_.facemap.find( key )->second.second;
      else
        return DGFBoundaryParameter::defaultValue();
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

    // return true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return dgf_.haveBndParameters;
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
    bool generateALUGrid( const ALUGridElementType eltype,
                          std::istream &file,
                          MPICommunicatorType communicator,
                          const std::string &filename );

    bool generateALU2dGrid( const ALUGridElementType eltype,
                            std::istream &file,
                            MPICommunicatorType communicator,
                            const std::string &filename );

    static Grid* callDirectly( const char* gridname,
                               const int rank,
                               const char *filename,
                               MPICommunicatorType communicator )
    {
  #if ALU3DGRID_PARALLEL
      // in parallel runs add rank to filename
      std :: stringstream tmps;
      tmps << filename << "." << rank;
      const std :: string &tmp = tmps.str();

      // if file exits then use it
      if( fileExists( tmp.c_str() ) )
        return new Grid( tmp.c_str(), communicator );
  #endif
      // for rank 0 we also check the normal file name
      if( rank == 0 )
      {
        if( fileExists( filename ) )
          return new Grid( filename , communicator );

        // only throw this exception on rank 0 because
        // for the other ranks we can still create empty grids
        DUNE_THROW( GridError, "Unable to create " << gridname << " from '"
                                                   << filename << "'." );
      }
      else
        dwarn << "WARNING:  P[" << rank << "]: Creating empty grid!" << std::endl;

      // return empty grid on all other processes
      return new Grid( communicator );
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

  // note: template parameter dimw is only added to avoid ALUSimplexGrid deprecation warning
  template < int dimw >
  struct DGFGridFactory< ALUSimplexGrid<3,dimw> > :
    public DGFBaseFactory< ALUSimplexGrid<3,dimw> >
  {
    typedef ALUSimplexGrid<3,dimw> DGFGridType;
    typedef DGFBaseFactory< DGFGridType > BaseType;
    typedef typename BaseType :: MPICommunicatorType MPICommunicatorType;
  protected:
    using BaseType :: grid_;
    using BaseType :: callDirectly;
  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< DGFGridType >( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< DGFGridType >( comm )
    {
      std::ifstream input( filename.c_str() );

      bool fileFound = input.is_open() ;
      if( fileFound )
        fileFound = generate( input, comm, filename );

      if( ! fileFound )
        grid_ = callDirectly( "ALUSimplexGrid< 3 , 3 >", this->rank( comm ), filename.c_str(), comm );
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  // note: template parameter dimw is only added to avoid ALUCubeGrid deprecation warning
  template < int dimw >
  struct DGFGridFactory< ALUCubeGrid<3, dimw> > :
    public DGFBaseFactory< ALUCubeGrid<3, dimw> >
  {
    typedef ALUCubeGrid<3,dimw> DGFGridType;
    typedef DGFBaseFactory< DGFGridType > BaseType;
    typedef typename BaseType :: MPICommunicatorType MPICommunicatorType;
  protected:
    using BaseType :: grid_;
    using BaseType :: callDirectly;
  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : DGFBaseFactory< DGFGridType >( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< DGFGridType >( comm )
    {
      std::ifstream input( filename.c_str() );
      bool fileFound = input.is_open() ;
      if( fileFound )
        fileFound = generate( input, comm, filename );

      if( ! fileFound )
        grid_ = callDirectly( "ALUCubeGrid< 3 , 3 >", this->rank( comm ), filename.c_str(), comm );
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  template < ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
  struct DGFGridFactory< ALUGrid<3,3, eltype, refinementtype, Comm > > :
    public DGFBaseFactory< ALUGrid<3,3, eltype, refinementtype, Comm > >
  {
    typedef ALUGrid<3,3, eltype, refinementtype, Comm > DGFGridType;
    typedef DGFBaseFactory< DGFGridType > BaseType;
    typedef typename BaseType :: MPICommunicatorType MPICommunicatorType;
  protected:
    using BaseType :: grid_;
    using BaseType :: callDirectly;
  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : BaseType( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : BaseType( comm )
    {
      std::ifstream input( filename.c_str() );
      bool fileFound = input.is_open() ;
      if( fileFound )
        fileFound = generate( input, comm, filename );

      if( ! fileFound )
        grid_ = callDirectly( "ALUGrid< 3, 3, eltype, ref, comm >", this->rank( comm ), filename.c_str(), comm );
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
  };

  template <int dimw>
  struct DGFGridFactory< ALUConformGrid<2,dimw> > :
    public DGFBaseFactory< ALUConformGrid<2,dimw> >
  {
    typedef DGFBaseFactory< ALUConformGrid<2,dimw> > Base;
    typedef typename Base:: MPICommunicatorType MPICommunicatorType;
    typedef typename Base::Grid Grid;
    const static int dimension = Grid::dimension;
    typedef typename Base::GridFactory GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : Base()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUConformGrid<2,dimw> >()
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( Base::fileExists( filename.c_str() ) )
          Base::grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }
  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
    using Base::grid_;
    using Base::factory_;
    using Base::dgf_;
  };

  template <int dimw>
  struct DGFGridFactory< ALUSimplexGrid<2,dimw> > :
    public DGFBaseFactory< ALUSimplexGrid<2,dimw> >
  {
    typedef DGFBaseFactory< ALUSimplexGrid<2,dimw> > Base;
    typedef typename Base::MPICommunicatorType MPICommunicatorType;
    typedef typename Base::Grid Grid;
    const static int dimension = Grid::dimension;
    typedef typename Base::GridFactory GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : Base()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW(DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUSimplexGrid<2,dimw> >()
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( Base::fileExists( filename.c_str() ) )
          Base::grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
    using Base::grid_;
    using Base::factory_;
    using Base::dgf_;
  };

  template <int dimw>
  struct DGFGridFactory< ALUCubeGrid<2,dimw> > :
    public DGFBaseFactory< ALUCubeGrid<2,dimw> >
  {
    typedef DGFBaseFactory< ALUCubeGrid<2,dimw> > Base;
    typedef typename Base::MPICommunicatorType MPICommunicatorType;
    typedef typename Base::Grid Grid;
    const static int dimension = Grid::dimension;
    typedef typename Base::GridFactory GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : Base()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW(DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : DGFBaseFactory< ALUCubeGrid<2,dimw> >()
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( Base::fileExists( filename.c_str() ) )
          Base::grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
    using Base::grid_;
    using Base::factory_;
    using Base::dgf_;
  };

  template < int dimw, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm >
  struct DGFGridFactory< ALUGrid<2, dimw, eltype, refinementtype, Comm > > :
    public DGFBaseFactory< ALUGrid< 2, dimw, eltype, refinementtype, Comm > >
  {
    typedef ALUGrid< 2, dimw, eltype, refinementtype, Comm > DGFGridType;
    typedef DGFBaseFactory< DGFGridType > BaseType;
    typedef typename BaseType :: MPICommunicatorType MPICommunicatorType;
    typedef typename BaseType::Grid Grid;
    const static int dimension = Grid::dimension;
    typedef typename BaseType::GridFactory GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
      : BaseType( comm )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW(DGFException, "Error resetting input stream." );
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator())
      : BaseType( comm )
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found." );
      if( !generate( input, comm, filename ) )
      {
        if( BaseType::fileExists( filename.c_str() ) )
          grid_ = new Grid( filename );
        else
          DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '" << filename << "'." );
      }
    }

  protected:
    bool generate( std::istream &file, MPICommunicatorType comm, const std::string &filename = "" );
    using BaseType::grid_;
    using BaseType::factory_;
    using BaseType::dgf_;
  };

}

#include "dgfalu.cc"

#endif // #if defined ENABLE_ALUGRID

#endif // #ifndef DUNE_DGFPARSERALU_HH
