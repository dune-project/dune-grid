// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
#define DUNE_ALU_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/utility/persistentcontainervector.hh>

namespace Dune
{

  // ALUGridPersistentContainer
  // --------------------------

  template< class G, class T >
  class ALUGridPersistentContainer
    : public PersistentContainerVector< G, typename G::HierarchicIndexSet, std::vector< T > >
  {
    typedef PersistentContainerVector< G, typename G::HierarchicIndexSet, std::vector< T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    ALUGridPersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid.hierarchicIndexSet(), codim, value )
    {}
  };


  // PersistentContainer for ALUGrid
  // -------------------------------

  template< int dim, int dimworld, class T >
  class PersistentContainer< ALUConformGrid< dim, dimworld >, T >
    : public ALUGridPersistentContainer< ALUConformGrid< dim, dimworld >, T >
  {
    typedef ALUGridPersistentContainer< ALUConformGrid< dim, dimworld >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, class T >
  class PersistentContainer< ALUCubeGrid< dim, dimworld >, T >
    : public ALUGridPersistentContainer< ALUCubeGrid< dim, dimworld >, T >
  {
    typedef ALUGridPersistentContainer< ALUCubeGrid< dim, dimworld >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, class T >
  class PersistentContainer< ALUSimplexGrid< dim, dimworld >, T >
    : public ALUGridPersistentContainer< ALUSimplexGrid< dim, dimworld >, T >
  {
    typedef ALUGridPersistentContainer< ALUSimplexGrid< dim, dimworld >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm, class T >
  class PersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T >
    : public ALUGridPersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T >
  {
    typedef ALUGridPersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, ALU2DSPACE ElementType elType, class T >
  class PersistentContainer< ALU2dGrid< dim, dimworld, elType >, T >
    : public ALUGridPersistentContainer< ALU2dGrid< dim, dimworld, elType >, T >
  {
    typedef ALUGridPersistentContainer< ALU2dGrid< dim, dimworld, elType >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

  template< ALU3dGridElementType elType, class Comm, class T >
  class PersistentContainer< ALU3dGrid< elType, Comm >, T >
    : public ALUGridPersistentContainer< ALU3dGrid< elType, Comm >, T >
  {
    typedef ALUGridPersistentContainer< ALU3dGrid< elType, Comm >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

} // namespace Dune

#endif // #if HAVE_ALUGRID

#endif // #ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
