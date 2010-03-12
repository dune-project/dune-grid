// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2D_ALUGRID_HH
#define DUNE_ALU2D_ALUGRID_HH

// only include this code, if ENABLE_ALUGRID is defined
#ifdef ENABLE_ALUGRID

#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

  /** @copydoc ALUSimplexGrid
      \brief [<em> provides \ref Dune::Grid </em>]
      \brief grid with support for simplicial mesh in 2d.
      \ingroup ALUSimplexGrid
   */
  template<int dimw>
  class ALUSimplexGrid< 2, dimw >
    : public Dune::ALU2dGrid< 2, dimw >
  {
    typedef ALUSimplexGrid< 2, dimw > This;

    typedef Dune::ALU2dGrid< 2, dimw > BaseType;
    enum { dim      = 2 };
    enum { dimworld = dimw };

  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUSimplexGrid(const std::string macroName,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0 )
      : BaseType(macroName,1, bndProject, bndVector)
    {
      std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }

    //! \brief constructor for creating ALUSimplexGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUSimplexGrid(const std::string macroName,
                   std::istream& macroFile,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0 )
      : BaseType("",1, bndProject, bndVector,&macroFile)
    {
      std::cout << "\nCreated serial ALUSimplexGrid<"<<dim<<","<<dimworld;
      if( macroName == "" )
        std::cout <<">. \n\n";
      else
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }

    //! constructor creating empty grid
    ALUSimplexGrid( ) : BaseType(1)
    {
      std::cout << "\nCreated empty ALUSimplexGrid<"<<dim<<","<<dimworld <<">. \n\n";
    }

    //! return name of the grid
    static inline std::string name () DUNE_DEPRECATED { return "ALUSimplexGrid"; }

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

    template< template< int, int > class, int >
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
   */
  template <int dim, int dimworld>
  class ALUConformGrid {};

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



  /** @copydoc Dune::ALUConformGrid
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d.
     \ingroup ALUConformGrid
   */
  template<int dimw>
  class ALUConformGrid< 2, dimw >
    : public Dune::ALU2dGrid< 2, dimw >
  {
    typedef ALUConformGrid< 2, dimw > This;

    typedef Dune::ALU2dGrid<2,dimw> BaseType;
    enum { dim      = 2 };
    enum { dimworld = dimw };
  public:
    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of boundary projection
    typedef typename BaseType :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

    //! \brief constructor for creating ALUConformGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUConformGrid(const std::string macroName,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0 )
      : BaseType(macroName, 0, bndProject, bndVector)
    {
      std::cout << "\nCreated serial ALUConformGrid<"<<dim<<","<<dimworld;
      std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }

    //! \brief constructor for creating ALUConformGrid from given macro grid file
    //! \param macroName filename for macro grid in ALUGrid triangle format
    ALUConformGrid(const std::string macroName,
                   std::istream& macroFile,
                   const DuneBoundaryProjectionType* bndProject  = 0,
                   const DuneBoundaryProjectionVector* bndVector = 0 )
      : BaseType("", 0, bndProject, bndVector, &macroFile)
    {
      std::cout << "\nCreated serial ALUConformGrid<"<<dim<<","<<dimworld;
      if( macroName == "" )
        std::cout <<">. \n\n";
      else
        std::cout <<"> from macro grid file '" << macroName << "'. \n\n";
    }

    //! constructor creating empty grid
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

    template< template< int, int > class, int >
    friend class ALU2dGridFactory;

    //! Copy constructor should not be used
    ALUConformGrid( const ALUConformGrid & g ) ; // : BaseType(g) {}

    //! assignment operator should not be used
    ALUConformGrid<dim,dimworld>&
    operator = (const ALUConformGrid& g);
  };

} //end  namespace Dune

#else
#error "Trying to use <dune/grid/alugrid.hh> without ALUGRID_CPPFLAGS."
#endif // #ifdef ENABLE_ALUGRID

#endif
