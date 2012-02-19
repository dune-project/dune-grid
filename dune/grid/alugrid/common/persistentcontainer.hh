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
  public:
    typedef ALUConformGrid< dim, dimworld > GridType;
  private:
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
  public:
    typedef ALUCubeGrid< dim, dimworld > GridType;
  private:
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
  public:
    typedef ALUSimplexGrid< dim, dimworld > GridType;
  private:
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

  template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm,
      class Data, class Allocator >
  class PersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, Data, Allocator >
    : public PersistentContainerVector< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >,
          typename ALUGrid< dim, dimworld, eltype, refinementtype, Comm >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
  public:
    typedef ALUGrid< dim, dimworld, eltype, refinementtype, Comm > GridType;
  private:
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

  template< int dim, int dimworld, ALU2DSPACE ElementType elType, class Data, class Allocator >
  class PersistentContainer< ALU2dGrid< dim, dimworld, elType >, Data, Allocator >
    : public PersistentContainerVector< ALU2dGrid< dim, dimworld, elType >,
          typename ALU2dGrid< dim, dimworld, elType >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
  public:
    typedef ALU2dGrid< dim, dimworld, elType >  GridType;
  private:
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

  template< ALU3dGridElementType elType, class Comm, class Data, class Allocator >
  class PersistentContainer< ALU3dGrid< elType, Comm >, Data, Allocator >
    : public PersistentContainerVector< ALU3dGrid< elType, Comm >,
          typename ALU3dGrid< elType, Comm >::HierarchicIndexSet,
          std::vector<Data,Allocator> >
  {
  public:
    typedef ALU3dGrid< elType, Comm >  GridType;
  private:
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  protected:
    using BaseType :: index_;
    using BaseType :: data_;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}

    //! this method is needed for the level communication
    //! of ALU3dGrid, see datahandle.hh
    const Data& getData ( const size_t idx ) const
    {
      assert( idx < data_.size() );
      return data_[ idx ];
    }
  };

} // end namespace Dune
#endif // ENABLE_ALU

#endif // end DUNE_ALU_PERSISTENTCONTAINER_HH
