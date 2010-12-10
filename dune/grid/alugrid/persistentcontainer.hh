// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
#define DUNE_ALU_PERSISTENTCONTAINER_HH

#ifdef ENABLE_ALUGRID
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/alugrid.hh>

namespace Dune
{
  // PersistentContainer for ALUGrid
  // -------------------------------

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUConformGrid< dim, dimworld >, Data, Allocator >
    : public PersistentContainerVector< ALUConformGrid< dim, dimworld >,
          typename ALUConformGrid< dim, dimworld >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
    typedef ALUConformGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUCubeGrid< dim, dimworld >, Data, Allocator >
    : public PersistentContainerVector< ALUCubeGrid< dim, dimworld >,
          typename ALUCubeGrid< dim, dimworld >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
    typedef ALUCubeGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUSimplexGrid< dim, dimworld >, Data, Allocator >
    : public PersistentContainerVector< ALUSimplexGrid< dim, dimworld >,
          typename ALUSimplexGrid< dim, dimworld >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
    typedef ALUSimplexGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };
} // end namespace Dune
#endif // ENABLE_ALU

#endif // end DUNE_ALU_PERSISTENTCONTAINER_HH
