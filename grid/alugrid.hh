// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// 3d version
#include <dune/grid/alu3dgrid.hh>

// 2d version
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune {

  template <int dim,int dimworld> class ALUCubeGrid;

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 3D grid with support for hexahedrons.
     @ingroup GridImplementations
     The ALUCubeGrid implements the Dune GridInterface for 3d hexahedral meshes.
     This grid can be locally adapted and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports hexahedrons and tetrahedrons.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)
   */
  template <>
  class ALUCubeGrid<3,3> :
    public Dune::ALU3dGrid<3,3,Dune::hexa> {
    typedef Dune::ALU3dGrid<3,3,Dune::hexa> BaseType;
    enum { dim      = 3 };
    enum { dimworld = 3 };
  public:
#if ALU3DGRID_PARALLEL
    //! constructor taking filename of macro grid and MPI_Comm
    ALUCubeGrid(const std::string macroName , MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm) {}
    //! constructor creating empty grid, empty string creates empty grid
    ALUCubeGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType("",mpiComm) {}
#else
    //! constructor taking filename of macro grid
    ALUCubeGrid(const std::string macroName , int mpicomm = 0 ) :
      BaseType(macroName) {}
    //! constructor creating empty grid
    ALUCubeGrid(int myrank = -1) :
      BaseType("",myrank) {}
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

    friend class Conversion< ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUCubeGrid<dimension,dimensionworld> , HasObjectStream > ;

  private:
    ALUCubeGrid(const ALUCubeGrid & g) : BaseType(g) {}
  };

  namespace Capabilities {

    template<int dim,int dimw>
    struct hasLeafIterator< Dune::ALUCubeGrid<dim, dimw > >
    {
      static const bool v = true;
    };

    template<int dim,int dimw, int cdim >
    struct hasEntity<Dune::ALUCubeGrid<dim, dimw>, cdim >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct isParallel<const ALUCubeGrid<dim, dimw> > {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct isLevelwiseConforming< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct hasHangingNodes< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct hasBackupRestoreFacilities< ALUCubeGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities


  template <int dim,int dimworld> class ALUSimplexGrid;

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     The ALUSimplexGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports hexahedrons and tetrahedrons.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)
   */
  template <>
  class ALUSimplexGrid<3,3> :
    public Dune::ALU3dGrid<3,3,Dune::tetra> {
    typedef Dune::ALU3dGrid<3,3,Dune::tetra> BaseType;
    enum { dim      = 3 };
    enum { dimworld = 3 };
  public:
#if ALU3DGRID_PARALLEL
    //! constructor taking filename of macro grid and MPI_Comm
    ALUSimplexGrid(const std::string macroName, MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm) {}
    //! constructor creating empty grid, empty string creates empty grid
    ALUSimplexGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType("",mpiComm) {}
#else
    //! constructor taking filename of macro grid
    ALUSimplexGrid(const std::string macroName , int mpicomm = 0) :
      BaseType(macroName) {}
    //! constructor creating empty grid
    ALUSimplexGrid(int myrank = -1) :
      BaseType("",myrank) {}
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

    friend class Conversion< ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;
    friend class Conversion< const ALUSimplexGrid<dimension,dimensionworld> , HasObjectStream > ;

  private:
    ALUSimplexGrid(const ALUSimplexGrid & g) : BaseType(g) {}
  };

  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief grid with support for simplicial mesh in 2d and 3d.
     @ingroup GridImplementations
     The ALUSimplexGrid implements the Dune GridInterface for 2d triangular and
     3d tetrahedral meshes.
     This grid can be locally adapted and used in parallel
     computations using dynamcic load balancing.

     @note
     Adaptive parallel grid supporting dynamic load balancing, written
     mainly by Bernard Schupp. This grid supports hexahedrons and tetrahedrons.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

     Two tools are available for partitioning :
     \li Metis ( version 4.0 and higher, see http://www-users.cs.umn.edu/~karypis/metis/metis/ )
     \li Party Lib ( version 1.1 and higher, see http://wwwcs.upb.de/fachbereich/AG/monien/RESEARCH/PART/party.html)
   */
  template <>
  class ALUSimplexGrid<2,2> :
    public Dune::ALU2dGrid<2,2> {
    typedef Dune::ALU2dGrid<2,2> BaseType;
    enum { dim      = 2 };
    enum { dimworld = 2 };
  public:
    ALUSimplexGrid(const std::string macroName
#if ALU3DGRID_PARALLEL
                   , MPI_Comm mpiComm = MPI_COMM_WORLD
#endif
                   ) :
      BaseType(macroName) {}
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
    typedef BaseType::HierarchicIteratorImp HierarchicIteratorImp;
  };

  namespace Capabilities {

    template<int dim,int dimw>
    struct hasLeafIterator< Dune::ALUSimplexGrid<dim, dimw > >
    {
      static const bool v = true;
    };

    template<int dim,int dimw, int cdim >
    struct hasEntity<Dune::ALUSimplexGrid<dim, dimw>, cdim >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct isParallel<const ALUSimplexGrid<dim, dimw> > {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct isLevelwiseConforming< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct hasHangingNodes< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    template<int dim,int dimw>
    struct hasBackupRestoreFacilities< ALUSimplexGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

} //end  namespace Dune
#endif
