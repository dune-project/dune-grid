// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2D_ALUGRID_HH
#define DUNE_ALU2D_ALUGRID_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/common/capabilities.hh>
#include <dune/grid/alugrid/2d/capabilities.hh>
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

  /** @copydoc ALUCubeGrid
      \brief [<em> provides \ref Dune::Grid </em>]
      \brief grid with support for cube mesh in 2d.
      \ingroup ALUCubeGrid
      \deprecated Use ALUGrid< 2, dimw, cube, nonconforming > instead.
   */
  template<int dimw>
  class ALUCubeGrid< 2, dimw >
    : public Dune::ALU2dGrid< 2, dimw, ALU2DSPACE quadrilateral >
  {
    typedef ALUCubeGrid< 2, dimw > This;

    typedef Dune::ALU2dGrid< 2, dimw, ALU2DSPACE quadrilateral > BaseType;
    enum { dim      = 2 };
    enum { dimworld = dimw };

  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, cube, nonconforming > instead.
    DUNE_DEPRECATED
    ALUCubeGrid(const std::string macroName,
                const DuneBoundaryProjectionType* bndProject  = 0,
                const DuneBoundaryProjectionVector* bndVector = 0,
                const bool verbose = true )
      : BaseType(macroName,1, bndProject, bndVector)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUCubeGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param macroFile  Stream to read macro grid file contents from.
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, cube, nonconforming > instead.
    DUNE_DEPRECATED
    ALUCubeGrid(const std::string macroName,
                std::istream& macroFile,
                const DuneBoundaryProjectionType* bndProject  = 0,
                const DuneBoundaryProjectionVector* bndVector = 0,
                const bool verbose = true )
      : BaseType("",1, bndProject, bndVector, &macroFile)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUCubeGrid<"<<dim<<","<<dimworld;
        if( macroName == "" )
          std::cout <<">. \n\n";
        else
          std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! constructor creating empty grid
    //! \deprecated Use ALUGrid< 2, dimw, cube, nonconforming > instead.
    DUNE_DEPRECATED
    ALUCubeGrid( ) : BaseType(1)
    {
      std::cout << "\nCreated empty ALUCubeGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 1 };
    typedef typename BaseType::ctype ctype;
    typedef typename BaseType::GridFamily GridFamily;
    typedef typename GridFamily::Traits Traits;
    typedef typename BaseType::LocalIdSetImp LocalIdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef typename BaseType::LeafIteratorImp LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;
    typedef typename BaseType::HierarchicIteratorImp HierarchicIteratorImp;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef typename Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
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
    friend class Conversion< ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUCubeGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUCubeGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    template< class >
    friend class ALU2dGridFactory;

    //! Copy constructor should not be used
    ALUCubeGrid( const ALUCubeGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUCubeGrid<dim,dimworld>&
    operator = (const ALUCubeGrid& g);
  };

  /** @copydoc ALUSimplexGrid
      \brief [<em> provides \ref Dune::Grid </em>]
      \brief grid with support for simplicial mesh in 2d.
      \ingroup ALUSimplexGrid
      \deprecated Use ALUGrid< 2, dimw, simplex, nonconforming > instead.
   */
  template<int dimw>
  class ALUSimplexGrid< 2, dimw >
    : public Dune::ALU2dGrid< 2, dimw, ALU2DSPACE triangle >
  {
    typedef ALUSimplexGrid< 2, dimw > This;

    typedef Dune::ALU2dGrid< 2, dimw, ALU2DSPACE triangle > BaseType;
    enum { dim      = 2 };
    enum { dimworld = dimw };

  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, simplex, nonconforming > instead.
    DUNE_DEPRECATED
    ALUSimplexGrid(const std::string macroName,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0,
                   const bool verbose = true )
      : BaseType(macroName,1, bndProject, bndVector)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param macroFile  Stream to read macro grid file contents from.
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, simplex, nonconforming > instead.
    DUNE_DEPRECATED
    ALUSimplexGrid(const std::string macroName,
                   std::istream& macroFile,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0,
                   const bool verbose = true )
      : BaseType("",1, bndProject, bndVector, &macroFile)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
        if( macroName == "" )
          std::cout <<">. \n\n";
        else
          std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! constructor creating empty grid
    //! \deprecated Use ALUGrid< 2, dimw, simplex, nonconforming > instead.
    DUNE_DEPRECATED
    ALUSimplexGrid( ) : BaseType(1)
    {
      std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 1 };
    typedef typename BaseType::ctype ctype;
    typedef typename BaseType::GridFamily GridFamily;
    typedef typename GridFamily::Traits Traits;
    typedef typename BaseType::LocalIdSetImp LocalIdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef typename BaseType::LeafIteratorImp LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;
    typedef typename BaseType::HierarchicIteratorImp HierarchicIteratorImp;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef typename Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
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
    friend class Conversion< ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUSimplexGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUSimplexGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    template< class >
    friend class ALU2dGridFactory;

    //! Copy constructor should not be used
    ALUSimplexGrid( const ALUSimplexGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUSimplexGrid<dim,dimworld>&
    operator = (const ALUSimplexGrid& g);
  };

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     @ingroup ALUConformGrid
     The ALUConformGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted (conforming) and used in parallel
     computations using dynamic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports triangular/tetrahedral elements -
     a 3d cube
     grid is also available via the grid implementation ALUCubeGrid or ALUSimplexGrid.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     \li Available Implementations
            - Dune::ALUConformGrid<2,dimw>

     For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .

     \note The unspecialized version of ALUConformGrid is defined for
           documentation purposes only, only the specializations can actually be
           used.
     \deprecated Use ALUGrid instead.
   */
  template <int dim, int dimworld>
  class ALUConformGrid {
    dune_static_assert((AlwaysFalse<integral_constant<int,dim> >::value),
                       "The unspecialized version of ALUConformGrid is defined "
                       "for documentation purposes only, only the "
                       "specializations can actually be used.");
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

    /** \brief ALUConformGrid has a conforming leaf grid
       \ingroup ALUConformGrid
     */
    template<int dim,int dimw>
    struct isLeafwiseConforming< ALUConformGrid<dim,dimw> >
    {
      static const bool v = true;
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



  /** @copydoc Dune::ALUConformGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d.
     \ingroup ALUConformGrid
     \deprecated Use ALUGrid< 2, dimw, simplex, conforming > instead.
   */
  template<int dimw>
  class ALUConformGrid< 2, dimw >
    : public Dune::ALU2dGrid< 2, dimw, ALU2DSPACE triangle >
  {
    typedef ALUConformGrid< 2, dimw > This;

    typedef Dune::ALU2dGrid<2,dimw, ALU2DSPACE triangle> BaseType;
    enum { dim      = 2 };
    enum { dimworld = dimw };
  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUConformGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, simplex, conforming > instead.
    DUNE_DEPRECATED
    ALUConformGrid(const std::string macroName,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0,
                   const bool verbose = true)
      : BaseType(macroName, 0, bndProject, bndVector)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUConformGrid<"<<dim<<","<<dimworld;
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! \brief constructor for creating ALUConformGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param macroFile  Stream to read macro grid file contents from.
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    //! \deprecated Use ALUGrid< 2, dimw, simplex, conforming > instead.
    DUNE_DEPRECATED
    ALUConformGrid(const std::string macroName,
                   std::istream& macroFile,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0,
                   const bool verbose = true )
      : BaseType("", 0, bndProject, bndVector, &macroFile)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial ALUConformGrid<"<<dim<<","<<dimworld;
        if( macroName == "" )
          std::cout <<">. \n\n";
        else
          std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
      }
    }

    //! constructor creating empty grid
    //! \deprecated Use ALUGrid< 2, dimw, simplex, conforming > instead.
    DUNE_DEPRECATED
    ALUConformGrid( ) : BaseType(0)
    {
      std::cout << "\nCreated empty ALUConformGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }

    //! return name of the grid
    static inline std::string name () { return "ALUConformGrid"; }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 2 };
    typedef typename BaseType::ctype ctype;
    typedef typename BaseType::GridFamily GridFamily;
    typedef typename GridFamily::Traits Traits;
    typedef typename BaseType::LocalIdSetImp LocalIdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef typename BaseType::LeafIteratorImp LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;
    typedef typename BaseType::HierarchicIteratorImp HierarchicIteratorImp;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef typename Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
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
    friend class Conversion< ALUConformGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUConformGrid<dimension,dimensionworld> , HasObjectStream > ;

    friend class Conversion< ALUConformGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;
    friend class Conversion< const ALUConformGrid<dimension,dimensionworld> , HasHierarchicIndexSet > ;

    template< class >
    friend class ALU2dGridFactory;

    //! Copy constructor should not be used
    ALUConformGrid( const ALUConformGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUConformGrid<dim,dimworld>&
    operator = (const ALUConformGrid& g);
  };

  /*-
     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     \li Available Implementations
          - quadrilateral and hexahedral elements only nonconforming refinement
            - Dune::ALUGrid< 2, 2, cube, nonconforming >
            - Dune::ALUGrid< 2, 3, cube, nonconforming >
            - Dune::ALUGrid< 3, 3, cube, nonconforming >
          - simplicial elements and nonconforming refinement
            - Dune::ALUGrid< 2, 2, simplex, nonconforming >
            - Dune::ALUGrid< 2, 3, simplex, nonconforming >
            - Dune::ALUGrid< 3, 3, simplex, nonconforming >
          - simplicial elements and bisection refinement
            - Dune::ALUGrid< 2, 2, simplex, conforming >
            - Dune::ALUGrid< 2, 3, simplex, conforming >
            - Dune::ALUGrid< 3, 3, simplex, conforming > (work in progress)

     \note template parameter Comm defaults to MPI_Comm, if MPI is available, No_Comm  otherwise.
   */
  template<int dimw, ALUGridElementType elType, ALUGridRefinementType refinementType, class Comm >
  class ALUGrid< 2, dimw, elType, refinementType, Comm >
    : public ALUGridBaseGrid < 2, dimw, elType, Comm > :: BaseGrid
  {
    typedef ALUGrid< 2, dimw, elType, refinementType, Comm > This;
    typedef typename ALUGridBaseGrid < 2, dimw, elType, Comm > :: BaseGrid BaseType;

    enum { dim      = 2 };
    enum { dimworld = dimw };

  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    ALUGrid(const std::string macroName,
            const DuneBoundaryProjectionType* bndProject  = 0,
            const DuneBoundaryProjectionVector* bndVector = 0,
            const bool verbose = true )
      : BaseType(macroName, hangingNodes(), bndProject, bndVector)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial " << name() << nameSuffix()
                  << " from macro grid file '" << macroName << "'." << std::endl << std::endl;
      }
    }

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    //! \param macroFile  Stream to read macro grid file contents from.
    //! \param bndProject global boundary projection pointer
    //! \param bndVector  pointer to vector holding boundary projection for
    //!                   each boundary segment.  ALUGrid takes ownership of
    //!                   this pointer and will delete it in the desctructor
    //! \param verbose    Whether to write a notice about grid creation to
    //!                   stdout.
    ALUGrid(const std::string macroName,
            std::istream& macroFile,
            const DuneBoundaryProjectionType* bndProject  = 0,
            const DuneBoundaryProjectionVector* bndVector = 0,
            const bool verbose = true )
      : BaseType("", hangingNodes(), bndProject, bndVector, &macroFile)
    {
      if( verbose )
      {
        std::cout << "\nCreated serial " << name() << nameSuffix();
        if( macroName != "" )
          std::cout <<" from macro grid file '" << macroName;
        std::cout << "." << std::endl << std::endl;
      }
    }

    static std::string name () { return std::string("ALUGrid"); }

    //! constructor creating empty grid
    ALUGrid( ) : BaseType( hangingNodes() )
    {
      std::cout << "\nCreated serial " << name() << nameSuffix() << "." << std::endl << std::endl;
    }

    enum {dimension=BaseType::dimension,dimensionworld=BaseType::dimensionworld};
    enum { refineStepsForHalf = 1 };
    typedef typename BaseType::ctype ctype;
    typedef typename BaseType::GridFamily GridFamily;
    typedef typename GridFamily::Traits Traits;
    typedef typename BaseType::LocalIdSetImp LocalIdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;
    typedef typename BaseType::LeafIteratorImp LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;
    typedef typename BaseType::HierarchicIteratorImp HierarchicIteratorImp;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Dune::GridView< DefaultLevelGridViewTraits< const This, pitype > >
      LevelGridView;
      typedef typename Dune::GridView< DefaultLeafGridViewTraits< const This, pitype > >
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
    static std::string nameSuffix()
    {
      std::string elt ( elType == cube ? "cube," : "simplex," );
      std::string ref ( refinementType == nonconforming ? "nonconforming>" : "conforming>" );
      std::stringstream suffix;
      suffix << "<"<<dim<<","<<dimworld<<"," << elt << ref;
      return suffix.str();
    }

    // returns number of hanging nodes allowed (0 or 1)
    int hangingNodes() const
    {
      return ((elType == simplex) && (refinementType == conforming)) ? 0 : 1;
    }

    friend class Conversion< This, HasObjectStream > ;
    friend class Conversion< const This, HasObjectStream > ;

    friend class Conversion< This, HasHierarchicIndexSet > ;
    friend class Conversion< const This, HasHierarchicIndexSet > ;

    template< class >
    friend class ALU2dGridFactory;

    //! Copy constructor should not be used
    ALUGrid( const ALUGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    This& operator = (const ALUGrid& g);
  };

} //end  namespace Dune

#else
#error "Trying to use <dune/grid/alugrid.hh> without ALUGRID_CPPFLAGS."
#endif // #ifdef ENABLE_ALUGRID

#endif
