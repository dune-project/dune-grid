// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// 3d version
#include <dune/grid/alu3dgrid.hh>

// 2d version
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune {

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
  template <int dim,int dimworld> class ALUCubeGrid;

  template <>
  class ALUCubeGrid<3,3> :
    public Dune::ALU3dGrid<3,3,Dune::hexa> {
    typedef Dune::ALU3dGrid<3,3,Dune::hexa> BaseType;
  public:
#ifdef _ALU3DGRID_PARALLEL_
    ALUCubeGrid(const std::string macroName , MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm) {}
    ALUCubeGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(mpiComm) {}
#else
    ALUCubeGrid(const std::string macroName ) :
      BaseType(macroName) {}
    ALUCubeGrid(int myrank = -1) :
      BaseType(myrank) {}
#endif
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
  template <int dim,int dimworld> class ALUSimplexGrid;

  template <>
  class ALUSimplexGrid<3,3> :
    public Dune::ALU3dGrid<3,3,Dune::tetra> {
    typedef Dune::ALU3dGrid<3,3,Dune::tetra> BaseType;
  public:
#ifdef _ALU3DGRID_PARALLEL_
    ALUSimplexGrid(const std::string macroName, MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(macroName,mpiComm) {}
    ALUSimplexGrid(MPI_Comm mpiComm = MPI_COMM_WORLD) :
      BaseType(mpiComm) {}
#else
    ALUSimplexGrid(const std::string macroName ) :
      BaseType(macroName) {}
    ALUSimplexGrid(int myrank = -1) :
      BaseType(myrank) {}
#endif
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
  public:
    ALUSimplexGrid(const std::string macroName
#ifdef _ALU3DGRID_PARALLEL_
                   , MPI_Comm mpiComm = MPI_COMM_WORLD
#endif
                   ) :
      BaseType(macroName) {}
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
