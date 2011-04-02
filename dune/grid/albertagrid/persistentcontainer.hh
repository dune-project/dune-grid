// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_PERSISTENTCONTAINER_HH
#define DUNE_ALBERTA_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // PersistentContainer for AlbertaGrid
  // -----------------------------------

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< AlbertaGrid< dim, dimworld >, Data, Allocator >
    : public PersistentContainerVector< AlbertaGrid< dim, dimworld >, typename AlbertaGrid< dim, dimworld >::HierarchicIndexSet, std::vector< Data, Allocator > >
  {
  public:
    typedef AlbertaGrid< dim, dimworld > GridType;
    private
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector< Data, Allocator > > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {}
  };

} // end namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
