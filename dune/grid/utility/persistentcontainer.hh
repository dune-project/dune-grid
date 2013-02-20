// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <map>

#include <dune/grid/utility/persistentcontainermap.hh>

namespace Dune
{

  /** \brief A class for storing data during an adaptation cycle.
   *
   * This container allows to store data which is to remain persistent
   * even during adaptation cycles. There is a default implementation based
   * on std::map but any grid implementation can provide a specialized implementation.
   *
   * The container expects that the class Data has a default constructor. This default
   * value is returned if data for an entity is requested which has not yet been added
   * to the container. Therefore no find as for std::map is required and provided.
   * The class furthermore provides a size method and iterators. The
   * guarantee being that the iteration is at least over all entries in the container
   * which are not set to default but the iteration could also include entries with the default
   * value. The size method returns the number of elements visited by an iteration so
   * that size() is larger or equal to the entries which are not set to the default value.
   * The method clear can be used to set all entries in the container to the default value.
   *
   * After the grid has changed update() or reserve() should to be called.
   * In contrast to the method reserve, the method update can compress the data to
   * reduce memory consumption, whereas reserve() can be implemented not to reduce any
   * memory even if the size of the grid has changed.
   *
   * Specialized implementations of this class are provided for some grid
   * implementations. Especially grids with an id type with a hash function could use
   * the std::unordered_map to store the data (the implementation Dune::PersistentContainerMap
   * can be used). For grids providing an id set suitable for storage in vector like
   * structures (id type is int and size method provided) Dune::PersistentContainerVector can be
   * used.
   */
  template< class G, class T >
  class PersistentContainer
    : public PersistentContainerMap< G, typename G::LocalIdSet, std::map< typename G::LocalIdSet::IdType, T > >
  {
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::map< typename G::LocalIdSet::IdType, T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, grid.localIdSet(), value )
    {}
  };

} // namespace Dune


#if 0

// the following implementation can be used for a grid providing a hash for the id type

#include <unordered_map>

namespace Dune
{

  template< G, class T >
  class PersistentContainer
    : public PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > >
  {
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value )
      : Base( grid, codim, grid.localIdSet(), value )
    {}
  };

} // namespace Dune

#endif // #if 0

namespace std
{

  template< class G, class T >
  inline void swap ( Dune::PersistentContainer< G, T > &a, Dune::PersistentContainer< G, T > &b )
  {
    a.swap( b );
  }

} // namespace std

#endif // #ifndef DUNE_PERSISTENTCONTAINER_HH
