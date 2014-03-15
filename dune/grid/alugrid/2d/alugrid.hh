// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2D_ALUGRID_HH
#define DUNE_ALU2D_ALUGRID_HH

// only include this code, if ENABLE_ALUGRID is defined
#if HAVE_ALUGRID || DOXYGEN

#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/common/capabilities.hh>
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

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

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelGridView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView LevelGridView;
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafGridView () const
    {
      typedef typename Partition< pitype >::LeafGridView LeafGridView;
      typedef typename LeafGridView::GridViewImp LeafGridViewImp;
      return LeafGridView( LeafGridViewImp( *this ) );
    }

    LevelGridView levelGridView ( int level ) const
    {
      typedef typename LevelGridView::GridViewImp LevelGridViewImp;
      return LevelGridView( LevelGridViewImp( *this, level ) );
    }

    LeafGridView leafGridView () const
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
#endif // #if HAVE_ALUGRID || DOXYGEN

#endif
