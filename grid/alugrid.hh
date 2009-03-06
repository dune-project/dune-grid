// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

// 3d version
#include "alugrid/3d/indexsets.hh"
#include "alugrid/3d/iterator.hh"
#include "alugrid/3d/entity.hh"
#include "alugrid/3d/geometry.hh"
#include "alugrid/3d/grid.hh"

// 2d version
#include <dune/grid/alugrid/2d/grid.hh>
/** @file
    @author Robert Kloefkorn
    @brief Provides base classes for ALUGrid
 **/

namespace Dune {


  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 3D grid with support for hexahedrons.
     @ingroup GridImplementations
     @ingroup ALUCubeGrid
     The ALUCubeGrid implements the Dune GridInterface for 3d hexahedral meshes.
     This grid can be locally adapted (non-conforming) and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports hexahedrons - a 2d/3d simplex
     grid is also available via the grid implementation ALUSimplexGrid or ALUConformGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)

     \li Available Implementations
          - Dune::ALUCubeGrid<3,3>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template< int dim, int dimworld >
  class ALUCubeGrid;

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
#if ALU3DGRID_PARALLEL
    //! \brief constructor for creating ALUCubeGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid hexa format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUCubeGrid(const std::string macroName , MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated parallel ALUCubeGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }
    //! \brief constructor creating empty grid
    ALUCubeGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType("",mpiComm)
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
    ALUCubeGrid(const std::string macroName , int mpiComm = 0 ) :
      BaseType(macroName)
    {
      std::cout << "\nCreated serial ALUCubeGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }
    //! constructor creating empty grid
    ALUCubeGrid(int myrank = -1) :
      BaseType("",myrank)
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
    ALUCubeGrid( const ALUCubeGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUCubeGrid<dim,dimworld>&
    operator = (const ALUCubeGrid& g);
  };

  namespace Capabilities {
    /** \struct isLeafwiseConforming
       \ingroup ALUCubeGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUCubeGrid
     */

    /** \brief ALUCubeGrid has entities for all codimension
       \ingroup ALUCubeGrid
     */
    template<int dim,int dimw, int cdim >
    struct hasEntity<Dune::ALUCubeGrid<dim, dimw>, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid is parallel
       \ingroup ALUCubeGrid
     */
    template<int dim,int dimw>
    struct isParallel<const ALUCubeGrid<dim, dimw> > {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid has conforming level grids
       \ingroup ALUCubeGrid
     */
    template<int dim,int dimw>
    struct isLevelwiseConforming< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid has supports hanging nodes
       \ingroup ALUCubeGrid
     */
    template<int dim,int dimw>
    struct hasHangingNodes< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief ALUCubeGrid has backup and restore facilities
       \ingroup ALUCubeGrid
     */
    template<int dim,int dimw>
    struct hasBackupRestoreFacilities< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities


  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUSimplexGrid
     The ALUSimplexGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted (non-conforming) and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports triangular/tetrahedral elements -
     a 3d cube
     grid is also available via the grid implementation ALUCubeGrid or ALUConformGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)

     \li Available Implementations
            - Dune::ALUSimplexGrid<3,3>
            - Dune::ALUSimplexGrid<2,2>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template< int dim, int dimworld >
  class ALUSimplexGrid;

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
#if ALU3DGRID_PARALLEL
    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid tetra format
    //! \param mpiComm MPI Communicator (when HAVE_MPI == 1 then mpiComm is of
    //!  type MPI_Comm and the default value is MPI_COMM_WORLD)
    ALUSimplexGrid(const std::string macroName, MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm)
    {
      if(this->comm().rank() == 0)
      {
        std::cout << "\nCreated parallel ALUSimplexGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }
    //! constructor creating empty grid, empty string creates empty grid
    ALUSimplexGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType("",mpiComm)
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
    ALUSimplexGrid(const std::string macroName , int mpicomm = 0) :
      BaseType(macroName)
    {
      std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }
    //! constructor creating empty grid
    ALUSimplexGrid(int myrank = -1) :
      BaseType("",myrank)
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
    ALUSimplexGrid( const ALUSimplexGrid & g ); //  : BaseType(g) {}

    //! assignment operator should not be used
    ALUSimplexGrid<dim,dimworld>&
    operator = (const ALUSimplexGrid& g);

  };

  /** @copydoc Dune::ALUSimplexGrid
      \brief [<em> provides \ref Dune::Grid </em>]
      \brief grid with support for simplicial mesh in 2d.
      \ingroup ALUSimplexGrid
   */
  template<>
  class ALUSimplexGrid< 2, 2 >
    : public Dune::ALU2dGrid< 2, 2 >
  {
    typedef ALUSimplexGrid< 2, 2 > This;

    typedef Dune::ALU2dGrid< 2, 2 > BaseType;
    enum { dim      = 2 };
    enum { dimworld = 2 };

  public:
    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUSimplexGrid(const std::string macroName )
      : BaseType(macroName,1)
    {
      std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }
    //! constructor creating empty grid
    ALUSimplexGrid( ) : BaseType(1)
    {
      std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }
    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 1 };
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
    typedef BaseType::HierarchicIteratorImp HierarchicIteratorImp;

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
    ALUSimplexGrid( const ALUSimplexGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUSimplexGrid<dim,dimworld>&
    operator = (const ALUSimplexGrid& g);
  };

  namespace Capabilities {
    /** \struct isLeafwiseConforming
       \ingroup ALUSimplexGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUSimplexGrid
     */

    /** \brief ALUSimplexGrid has entities for all codimension
       \ingroup ALUSimplexGrid
     */
    template<int dim,int dimw, int cdim >
    struct hasEntity<Dune::ALUSimplexGrid<dim, dimw>, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid is parallel
       \ingroup ALUSimplexGrid
     */
    template<int dim,int dimw>
    struct isParallel<const ALUSimplexGrid<dim, dimw> > {
      static const bool v = false;
    };

    /** \brief ALUSimplexGrid has conforming level grids
       \ingroup ALUSimplexGrid
     */
    template<int dim,int dimw>
    struct isLevelwiseConforming< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid has supports hanging nodes
       \ingroup ALUSimplexGrid
     */
    template<int dim,int dimw>
    struct hasHangingNodes< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief ALUSimplexGrid has backup and restore facilities
       \ingroup ALUSimplexGrid
     */
    template<int dim,int dimw>
    struct hasBackupRestoreFacilities< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUConformGrid
     The ALUConformGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted (conforming) and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports triangular/tetrahedral elements -
     a 3d cube
     grid is also available via the grid implementation ALUCubeGrid or ALUSimplexGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     \li Available Implementations
            - Dune::ALUConformGrid<2,2>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
   */
  template <int dim, int dimworld>
  class ALUConformGrid {};

  /** @copydoc Dune::ALUConformGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d.
     \ingroup ALUConformGrid
   */
  template<>
  class ALUConformGrid< 2, 2 >
    : public Dune::ALU2dGrid< 2, 2 >
  {
    typedef ALUConformGrid< 2, 2 > This;

    typedef Dune::ALU2dGrid<2,2> BaseType;
    enum { dim      = 2 };
    enum { dimworld = 2 };
  public:
    //! \brief constructor for creating ALUConformGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUConformGrid(const std::string macroName )
      : BaseType(macroName)
    {
      std::cout << "\nCreated serial ALUConformGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }
    //! constructor creating empty grid
    ALUConformGrid( ) : BaseType(0)
    {
      std::cout << "\nCreated empty ALUConformGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }
    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 2 };
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
    typedef BaseType::HierarchicIteratorImp HierarchicIteratorImp;

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
    friend class Conversion< ALUConformGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUConformGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUConformGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUConformGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    //! Copy constructor should not be used
    ALUConformGrid( const ALUConformGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUConformGrid<dim,dimworld>&
    operator = (const ALUConformGrid& g);
  };

  namespace Capabilities {
    /** \struct isLeafwiseConforming
       \ingroup ALUConformGrid
     */

    /** \struct IsUnstructured
       \ingroup ALUConformGrid
     */

    /** \brief ALUConformGrid has entities for all codimension
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw, int cdim >
    struct hasEntity<Dune::ALUConformGrid<dim, dimw>, cdim >
    {
      static const bool v = true;
    };

    /** \brief ALUConformGrid is parallel
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw>
    struct isParallel<const ALUConformGrid<dim, dimw> > {
      static const bool v = false;
    };

    /** \brief ALUConformGrid has non-conforming level grids
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw>
    struct isLevelwiseConforming< ALUConformGrid<dim,dimw> >
    {
      static const bool v = false;
    };

    /** \brief ALUConformGrid does not support hanging nodes
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw>
    struct hasHangingNodes< ALUConformGrid<dim,dimw> >
    {
      static const bool v = false;
    };

    /** \brief ALUConformGrid has backup and restore facilities
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw>
    struct hasBackupRestoreFacilities< ALUConformGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities


} //end  namespace Dune

#else
#error "Trying to use <dune/grid/alugrid.hh> without ALUGRID_CPPFLAGS."
#endif // #ifdef ENABLE_ALUGRID

#endif
