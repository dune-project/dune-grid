// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// 3d version
#include <dune/grid/alu3dgrid.hh>

// 2d version
#include <dune/grid/alu2dgrid.hh>

namespace Dune {

  // ALUGrid version supporting cubes
  template <int dim,int dimworld> class ALUCubeGrid;

  // ALUGrid version supporting simplices
  template <int dim,int dimworld> class ALUSimplexGrid;

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

  template <>
  class ALUSimplexGrid<2,2> :
    public Dune::ALU2dGrid<2,2> {
    typedef Dune::ALU2dGrid<2,2> BaseType;
  public:
    ALUSimplexGrid(const std::string macroName ) :
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

} //end  namespace Dune
#endif
