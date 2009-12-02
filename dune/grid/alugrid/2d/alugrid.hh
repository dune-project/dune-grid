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

    //! return name of the grid
    static inline std::string name () { return "ALUSimplexGrid"; }

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
            - Dune::ALUConformGrid<2,2>

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

    //! return name of the grid
    static inline std::string name () { return "ALUConformGrid"; }

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

} //end  namespace Dune

#else
#error "Trying to use <dune/grid/alugrid.hh> without ALUGRID_CPPFLAGS."
#endif // #ifdef ENABLE_ALUGRID

#endif
