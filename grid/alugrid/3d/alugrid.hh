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



  /** @copydoc Dune::ALUCubeGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 3D grid with support for hexahedrons.
     @ingroup ALUCubeGrid
   */
  template<>
  class ALUCubeGrid< 3, 3 >
    : public Dune::ALU3dGrid< 3, 3, Dune::hexa >
  {
    typedef ALUCubeGrid< 3, 3 > This;

    typedef Dune::ALU3dGrid<3,3,Dune::hexa> BaseType;
    enum { dim      = 3 };
    enum { dimworld = 3 };
  public:
    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

  #if ALU3DGRID_PARALLEL
    //! \brief constructor for creating ALUCubeGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid hexa format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUCubeGrid(const std::string macroName,
                MPI_Comm mpiComm = MPI_COMM_WORLD,
                const DuneBoundaryProjectionType* bndProject = 0) :
      BaseType(macroName,mpiComm,bndProject)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated parallel ALUCubeGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }
    //! \brief constructor creating empty grid
    ALUCubeGrid(MPI_Comm mpiComm = MPI_COMM_WORLD,
                const DuneBoundaryProjectionType* bndProject = 0) :
      BaseType("",mpiComm, bndProject)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated empty ALUCubeGrid<"<<dim<<","<<dimworld <<">. \n\n";
      }
    }
  #else
    //! \brief constructor for creating ALUCubeGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid hexa format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUCubeGrid(const std::string macroName,
                const DuneBoundaryProjectionType* bndProject = 0)
      : BaseType(macroName, bndProject)
    {
      std::cout << "\nCreated serial ALUCubeGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }

    //! constructor creating empty grid
    ALUCubeGrid() : BaseType("", (const DuneBoundaryProjectionType*) 0)
    {
      std::cout << "\nCreated empty ALUCubeGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }
  #endif
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

    //! Copy constructor should not be used
    ALUCubeGrid( const ALUCubeGrid & g ) ;   // : BaseType(g) {}

    //! assignment operator should not be used
    ALUCubeGrid<dim,dimworld>&
    operator = (const ALUCubeGrid& g);
  };



  /**  @copydoc Dune::ALUSimplexGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 3d.
     \ingroup ALUSimplexGrid
   */
  template<>
  class ALUSimplexGrid< 3, 3 >
    : public Dune::ALU3dGrid< 3, 3,Dune::tetra >
  {
    typedef ALUSimplexGrid< 3, 3 > This;

    typedef Dune::ALU3dGrid<3,3,Dune::tetra> BaseType;
    enum { dim      = 3 };
    enum { dimworld = 3 };
  public:
    //! type of boundary projection
    typedef BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

  #if ALU3DGRID_PARALLEL
    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid tetra format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUSimplexGrid(const std::string macroName,
                   MPI_Comm mpiComm = MPI_COMM_WORLD,
                   const DuneBoundaryProjectionType* bndProject = 0) :
      BaseType(macroName, mpiComm, bndProject)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated parallel ALUSimplexGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }
    //! constructor creating empty grid, empty string creates empty grid
    ALUSimplexGrid(MPI_Comm mpiComm = MPI_COMM_WORLD,
                   const DuneBoundaryProjectionType* bndProject = 0) :
      BaseType("", mpiComm, bndProject)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
      }
    }
  #else
    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid tetra format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUSimplexGrid(const std::string macroName,
                   const DuneBoundaryProjectionType* bndProject = 0) :
      BaseType(macroName, bndProject)
    {
      std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }
    //! constructor creating empty grid
    ALUSimplexGrid() :
      BaseType("", (const DuneBoundaryProjectionType*) 0)
    {
      std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }
  #endif
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

    //! Copy constructor should not be used
    ALUSimplexGrid( const ALUSimplexGrid & g );   //  : BaseType(g) {}

    //! assignment operator should not be used
    ALUSimplexGrid<dim,dimworld>&
    operator = (const ALUSimplexGrid& g);
  };

} //end  namespace Dune

#endif // #ifdef ENABLE_ALUGRID

#endif
