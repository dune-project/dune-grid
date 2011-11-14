// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_ALUGRID_HH
#define DUNE_ALU3DGRID_ALUGRID_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

// 3d version
#include <dune/grid/alugrid/3d/capabilities.hh>
#include <dune/grid/alugrid/3d/indexsets.hh>
#include <dune/grid/alugrid/3d/iterator.hh>
#include <dune/grid/alugrid/3d/entity.hh>
#include <dune/grid/alugrid/3d/geometry.hh>
#include <dune/grid/alugrid/3d/grid.hh>

/** @file
    @author Robert Kloefkorn
    @brief Provides base classes for ALUGrid
 **/

namespace Dune
{

  template< int dim, int dimworld >
  class ALUCubeGrid;

  template< int dim, int dimworld >
  class ALUSimplexGrid;

  static const char* ALUGridParallelSerial()
  {
#if ALU3DGRID_PARALLEL
    return "parallel";
#else
    return "serial";
#endif
  }

  /** @copydoc ALUCubeGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 3D grid with support for hexahedrons.
     @ingroup ALUCubeGrid
   */
  template<>
  class ALUCubeGrid< 3, 3 >
    : public Dune::ALU3dGrid< hexa >
  {
    typedef ALUCubeGrid< 3, 3 > This;
    typedef Dune::ALU3dGrid< hexa > BaseType;

    enum { dim      = 3 };
    enum { dimworld = 3 };

    typedef BaseType::MPICommunicatorType MPICommunicatorType;

  public:
    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUCubeGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid hexa format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verb       Whether to write a notice about grid creation to
    //!                   stdout.
    ALUCubeGrid(const std::string macroName,
                const MPICommunicatorType mpiComm = BaseType::defaultCommunicator(),
                const DuneBoundaryProjectionType* bndProject = 0,
                const DuneBoundaryProjectionVector* bndVector= 0,
                const bool verb = true ) :
      BaseType(macroName,mpiComm,bndProject, bndVector, nonconforming )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " ALUCubeGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor called from ALUGridFactory
    //! for creating ALUCubeGrid from given macro grid file
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of type MPI_Comm)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for each boundary segment
    //! \note ALUGrid takes ownership of this pointer and will delete it in the desctructor
    //! \param macroName filename from which ALUGrid is being generated
    //! \param verb      Whether to write a notice about grid creation to
    //!                  stdout.
    ALUCubeGrid(const MPICommunicatorType mpiComm,
                const DuneBoundaryProjectionType* bndProject ,
                const DuneBoundaryProjectionVector* bndVector,
                const std::string macroName,
                const bool verb = true ) :
      BaseType("", mpiComm, bndProject, bndVector, nonconforming )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " ALUCubeGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor creating empty grid
    ALUCubeGrid(const MPICommunicatorType mpiComm = BaseType::defaultCommunicator() ) :
      BaseType("", mpiComm,
               (const DuneBoundaryProjectionType *) 0,
               (const DuneBoundaryProjectionVector* ) 0,
               nonconforming )
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated empty ALUCubeGrid<"<<dim<<","<<dimworld <<">. \n\n";
      }
    }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    typedef BaseType::ctype ctype;
    typedef BaseType::GridFamily GridFamily;
    typedef GridFamily::Traits Traits;
    typedef BaseType::LocalIdSetImp LocalIdSetImp;
    typedef Traits :: GlobalIdSet GlobalIdSet;
    typedef Traits :: LocalIdSet LocalIdSet;
    typedef GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef BaseType::LeafIteratorImp LeafIteratorImp;
    typedef Traits::Codim<0>::LeafIterator LeafIteratorType;
    typedef Traits::Codim<0>::LeafIterator LeafIterator;

    // ALUGrid only typedefs
    typedef BaseType::HierarchicIteratorImp HierarchicIteratorImp;
    typedef BaseType::ObjectStreamType ObjectStreamType;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
      LeafGridView;
    };

    typedef Partition< All_Partition > :: LevelGridView LevelGridView;
    typedef Partition< All_Partition > :: LeafGridView LeafGridView;

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView LevelGridView;
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafView () const
    {
      typedef typename Partition< pitype >::LeafGridView LeafGridView;
      typedef typename LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

    LevelGridView levelView ( int level ) const
    {
      typedef LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    LeafGridView leafView () const
    {
      typedef LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

  private:
    friend class Conversion< ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUCubeGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUCubeGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    template< class > friend class ALU3dGridFactory;

    //! Copy constructor should not be used
    ALUCubeGrid( const ALUCubeGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUCubeGrid<dim,dimworld>&
    operator = (const ALUCubeGrid& g);
  };



  /**  @copydoc ALUSimplexGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 3d.
     \ingroup ALUSimplexGrid
   */
  template<>
  class ALUSimplexGrid< 3, 3 >
    : public Dune::ALU3dGrid< tetra >
  {
    typedef ALUSimplexGrid< 3, 3 > This;
    typedef Dune::ALU3dGrid< tetra > BaseType;

    enum { dim      = 3 };
    enum { dimworld = 3 };

    typedef BaseType::MPICommunicatorType MPICommunicatorType;

  public:
    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid tetra format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verb       Whether to write a notice about grid creation to
    //!                   stdout.
    ALUSimplexGrid(const std::string macroName,
                   const MPICommunicatorType mpiComm = BaseType::defaultCommunicator(),
                   const DuneBoundaryProjectionType* bndProject = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0,
                   const bool verb = true ) :
      BaseType(macroName, mpiComm, bndProject, bndVector, nonconforming )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " ALUSimplexGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor called from ALUGridFactory
    //! for creating ALUSimplexGrid from given macro grid file
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of type MPI_Comm)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for each boundary segment
    //!  \note ALUGrid takes ownership of this pointer and will delete it in the desctructor
    //! \param macroName filename from which ALUGrid is being generated
    //! \param verb       Whether to write a notice about grid creation to
    //!                   stdout.
    ALUSimplexGrid(const MPICommunicatorType mpiComm,
                   const DuneBoundaryProjectionType* bndProject ,
                   const DuneBoundaryProjectionVector* bndVector,
                   const std::string macroName,
                   const bool verb = true ) :
      BaseType("", mpiComm, bndProject, bndVector, nonconforming )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " ALUSimplexGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! constructor creating empty grid, empty string creates empty grid
    ALUSimplexGrid(const MPICommunicatorType mpiComm = BaseType::defaultCommunicator()) :
      BaseType("", mpiComm,
               (const DuneBoundaryProjectionType *) 0,
               (const DuneBoundaryProjectionVector* ) 0,
               nonconforming )
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
      }
    }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    typedef BaseType::ctype ctype;
    typedef BaseType::GridFamily GridFamily;
    typedef GridFamily::Traits Traits;
    typedef BaseType::LocalIdSetImp LocalIdSetImp;
    typedef Traits :: GlobalIdSet GlobalIdSet;
    typedef Traits :: LocalIdSet LocalIdSet;
    typedef GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef BaseType::LeafIteratorImp LeafIteratorImp;
    typedef Traits::Codim<0>::LeafIterator LeafIteratorType;
    typedef Traits::Codim<0>::LeafIterator LeafIterator;

    // ALUGrid only typedefs
    typedef BaseType::HierarchicIteratorImp HierarchicIteratorImp;
    typedef BaseType::ObjectStreamType ObjectStreamType;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
      LeafGridView;
    };

    typedef Partition< All_Partition > :: LevelGridView LevelGridView;
    typedef Partition< All_Partition > :: LeafGridView LeafGridView;

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView LevelGridView;
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafView () const
    {
      typedef typename Partition< pitype >::LeafGridView LeafGridView;
      typedef typename LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

    LevelGridView levelView ( int level ) const
    {
      typedef LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    LeafGridView leafView () const
    {
      typedef LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

  private:
    friend class Conversion< ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUSimplexGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUSimplexGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    template< class > friend class ALU3dGridFactory;

    //! Copy constructor should not be used
    ALUSimplexGrid( const ALUSimplexGrid & g ); //  : BaseType(g) {}

    //! assignment operator should not be used
    ALUSimplexGrid<dim,dimworld>&
    operator = (const ALUSimplexGrid& g);
  };

  /**  @copydoc ALUGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with 3d implementations of the grid interface
     \ingroup ALUGrid
   */
  template< ALUGridElementType elType, ALUGridRefinementType refineType >
  class ALUGrid< 3, 3, elType, refineType >
    : public ALUGridBaseGrid< elType > :: BaseGrid
  {
    typedef ALUGrid< 3, 3, elType, refineType > This;
    typedef typename ALUGridBaseGrid< elType > :: BaseGrid BaseType;

    enum { dim      = 3 };
    enum { dimworld = 3 };

    typedef typename BaseType::MPICommunicatorType MPICommunicatorType;

  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid tetra format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verb       Whether to write a notice about grid creation to
    //!                   stdout.
    ALUGrid(const std::string macroName,
            const MPICommunicatorType mpiComm = BaseType::defaultCommunicator(),
            const DuneBoundaryProjectionType* bndProject = 0,
            const DuneBoundaryProjectionVector* bndVector = 0,
            const bool verb = true ) :
      BaseType(macroName, mpiComm, bndProject, bndVector, refineType )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " " << name() << nameSuffix()
                  << " from macro grid file '" << macroName << "'. \n\n";
      }
    }

    static std::string name () { return std::string("ALUGrid"); }

    static std::string nameSuffix()
    {
      std::string elt ( elType == cube ? "cube," : "simplex," );
      std::string ref ( refineType == nonconforming ? "nonconforming>" : "conforming>" );
      std::stringstream suffix;
      suffix << "<"<<dim<<","<<dimworld<<"," << elt << ref;
      return suffix.str();
    }


    //! \brief constructor called from ALUGridFactory
    //! for creating ALUConformGrid from given macro grid file
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of type MPI_Comm)
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for each boundary segment
    //!  \note ALUGrid takes ownership of this pointer and will delete it in the desctructor
    //! \param macroName filename from which ALUGrid is being generated
    //! \param verb       Whether to write a notice about grid creation to
    //!                   stdout.
    ALUGrid(const MPICommunicatorType mpiComm,
            const DuneBoundaryProjectionType* bndProject ,
            const DuneBoundaryProjectionVector* bndVector,
            const std::string macroName,
            const bool verb = true ) :
      BaseType("", mpiComm, bndProject, bndVector, refineType )
    {
      const bool verbose = verb && this->comm().rank() == 0;
      if( verbose )
      {
        std::cout << "\nCreated " << ALUGridParallelSerial() << " " << name() << nameSuffix()
                  << " from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! constructor creating empty grid, empty string creates empty grid
    ALUGrid(const MPICommunicatorType mpiComm = BaseType::defaultCommunicator()) :
      BaseType("", mpiComm,
               (const DuneBoundaryProjectionType *) 0,
               (const DuneBoundaryProjectionVector* ) 0,
               refineType )
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated empty " << ALUGridParallelSerial() << " " << name() << nameSuffix() << "." << std::endl;
      }
    }

    enum { dimension=BaseType::dimension,  dimensionworld=BaseType::dimensionworld};
    typedef typename BaseType::ctype ctype;
    typedef typename BaseType::GridFamily GridFamily;
    typedef typename GridFamily::Traits Traits;
    typedef typename BaseType::LocalIdSetImp LocalIdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef typename BaseType::LeafIteratorImp LeafIteratorImp;
    typedef typename Traits:: template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits:: template Codim<0>::LeafIterator LeafIterator;

    // ALUGrid only typedefs
    typedef typename BaseType::HierarchicIteratorImp HierarchicIteratorImp;
    typedef typename BaseType::ObjectStreamType ObjectStreamType;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
      LeafGridView;
    };

    typedef typename Partition< All_Partition > :: LevelGridView LevelGridView;
    typedef typename Partition< All_Partition > :: LeafGridView LeafGridView;

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView LevelGridView;
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafView () const
    {
      typedef typename Partition< pitype >::LeafGridView LeafGridView;
      typedef typename LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

    LevelGridView levelView ( int level ) const
    {
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    LeafGridView leafView () const
    {
      typedef typename LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

  private:
    friend class Conversion< This , HasObjectStream > ;
    friend class Conversion< const This, HasObjectStream > ;

    friend class Conversion< This, HasHierarchicIndexSet > ;
    friend class Conversion< const This, HasHierarchicIndexSet > ;

    template< class > friend class ALU3dGridFactory;

    //! Copy constructor should not be used
    ALUGrid( const ALUGrid & g ); //  : BaseType(g) {}

    //! assignment operator should not be used
    This& operator = (const ALUGrid& g);
  };

} //end  namespace Dune

#endif // #ifdef ENABLE_ALUGRID

#undef alu_inline
#endif
