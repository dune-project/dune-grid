// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH
#define DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

#ifdef DOXYGEN

  //! Iterates over all entities of a GridView with the given codimension.
  /**
   * \relates GridView
   */
  template<typename GV, int codim, unsigned int partitions>
  IteratorRange<...> entities(const GV& gv, Codim<codim>, PartitionSet<partitions>);

  //! Iterates over all entities of a GridView with the given codimension.
  /**
   * \relates GridView
   */
  template<typename GV, int codim>
  IteratorRange<...> entities(const GV& gv, Codim<codim> cd);

  //! Iterates over all entities of a GridView with the given dimension.
  /**
   * \relates GridView
   */
  template<typename GV, int dim, unsigned int partitions>
  IteratorRange<...> entities(const GV& gv, Dim<dim>, PartitionSet<partitions>);

  //! Iterates over all entities of a GridView with the given dimension.
  /**
   * \relates GridView
   */
  template<typename GV, int dim>
  IteratorRange<...> entities(const GV& gv, Dim<dim>);

  //! Iterates over all elements (entities with codimension 0) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> elements(const GV& gv, PartitionSet<partitions>);

  //! Iterates over all elements (entities with codimension 0) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV>
  IteratorRange<...> elements(const GV& gv);

  //! Iterates over all facets (entities with codimension 1) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> facets(const GV& gv, PartitionSet<partitions>);

  //! Iterates over all facets (entities with codimension 1) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV>
  IteratorRange<...> facets(const GV& gv);

  //! Iterates over all edges (entities with dimension 1) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> edges(const GV& gv, PartitionSet<partitions>);

  //! Iterates over all edges (entities with dimension 1) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV>
  IteratorRange<...> edges(const GV& gv);

  //! Iterates over all vertices (entities with dimension 0) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> vertices(const GV& gv, PartitionSet<partitions>);

  //! Iterates over all vertices (entities with dimension 0) of a GridView.
  /**
   * \relates GridView
   */
  template<typename GV>
  IteratorRange<...> vertices(const GV& gv);

#else // DOXYGEN

  template<typename GV, int codim, unsigned int partitions>
  IteratorRange<
    typename GV::template Codim<codim>::template Partition<
      derive_partition_iterator_type<partitions>::value
      >::Iterator
    >
  entities(const GV& gv, Codim<codim>, PartitionSet<partitions>)
  {
    static_assert(0 <= codim && codim <= GV::dimension, "invalid codimension for given GridView");
    const PartitionIteratorType pit = derive_partition_iterator_type<partitions>::value;
    typedef IteratorRange<
      typename GV::template Codim<codim>::template Partition<pit>::Iterator
      > return_type;
    return return_type(gv.template begin<codim,pit>(),gv.template end<codim,pit>());
  }

  template<typename GV, int codim>
  auto entities(const GV& gv, Codim<codim> cd)
    -> decltype(entities(gv,cd,Partitions::all))
  {
    static_assert(0 <= codim && codim <= GV::dimension, "invalid codimension for given GridView");
    return entities(gv,cd,Partitions::all);
  }

  template<typename GV, int dim, unsigned int partitions>
  auto entities(const GV& gv, Dim<dim>, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<GV::dimension - dim>(),PartitionSet<partitions>()))
  {
    static_assert(0 <= dim && dim <= GV::dimension, "invalid dimension for given GridView");
    return entities(gv,Codim<GV::dimension - dim>(),PartitionSet<partitions>());
  }

  template<typename GV, int dim>
  auto entities(const GV& gv, Dim<dim>)
    -> decltype(entities(gv,Codim<GV::dimension - dim>()))
  {
    static_assert(0 <= dim && dim <= GV::dimension, "invalid dimension for given GridView");
    return entities(gv,Codim<GV::dimension - dim>());
  }

  template<typename GV, unsigned int partitions>
  auto elements(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<0>(),PartitionSet<partitions>()))
  {
    return entities(gv,Codim<0>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto elements(const GV& gv)
    -> decltype(entities(gv,Codim<0>()))
  {
    return entities(gv,Codim<0>());
  }

  template<typename GV, unsigned int partitions>
  auto facets(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<1>(),PartitionSet<partitions>()))
  {
    return entities(gv,Codim<1>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto facets(const GV& gv)
    -> decltype(entities(gv,Codim<1>()))
  {
    return entities(gv,Codim<1>());
  }

  template<typename GV, unsigned int partitions>
  auto edges(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Dim<1>(),PartitionSet<partitions>()))
  {
    return entities(gv,Dim<1>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto edges(const GV& gv)
    -> decltype(entities(gv,Dim<1>()))
  {
    return entities(gv,Dim<1>());
  }

  template<typename GV, unsigned int partitions>
  auto vertices(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Dim<0>(),PartitionSet<partitions>()))
  {
    return entities(gv,Dim<0>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto vertices(const GV& gv)
    -> decltype(entities(gv,Dim<0>()))
  {
    return entities(gv,Dim<0>());
  }

#endif // DOXYGEN

  //! Iterates over all intersections of an Entity with respect to the given GridView.
  /**
   * \relates GridView
   * \relates Entity
   */
  template<typename GV, typename Entity>
  IteratorRange<typename GV::IntersectionIterator> intersections(const GV& gv, const Entity& e)
  {
    return IteratorRange<typename GV::IntersectionIterator>(gv.ibegin(e),gv.iend(e));
  }

} // namespace Dune

#endif // DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH
