// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_VERTEXORDERFACTORY_HH
#define DUNE_GRID_UTILITY_VERTEXORDERFACTORY_HH

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/generalvertexorder.hh>

namespace Dune {

  //! Factory for GeneralVertexOrder objects using an IdSet
  /**
   * \tparam IdSet Type used to get the ids of the vertices.
   * \tparam Index Type of the indices provided by the vertex ordering
   *               object.  Must be integral, may be non-negative.
   *
   * \warning The Interface of the VertexOrder stuff is subject to change.  It
   *          is currently needed to use some global-valued finite elements
   *          from dune-localfunctions.
   *
   * \sa GeneralVertexOrder, reduceOrder()
   */
  template<class IdSet, class Index = std::size_t>
  class VertexOrderByIdFactory {
    const IdSet& idset;

  public:
    //! type of vertex order object may depend on the dimension of the element
    template<std::size_t dim>
    struct VertexOrder {
      //! type of vertex order object
      typedef GeneralVertexOrder<dim, Index> type;
    };

    //! construct a factory object
    /**
     * \tparam idset_ IdSet to use to extract the vertex ids.
     *
     * This factory object stores a reference to the IdSet object.  The
     * factory object's value will become singular when the stored reference
     * becomes invalid.  The only valid operation on a factory with singular
     * value is destruction, all other operations will result in undefined
     * behaviour.
     */
    VertexOrderByIdFactory(const IdSet &idset_) : idset(idset_) { }

    //! construct a vertex ordering object
    /**
     * \param e Grid element to create the vertex ordering object for.
     *
     * The returned object will remain valid even after the factory has become
     * singular or has been destroyed.
     */
    template<typename Element>
    typename VertexOrder<Element::mydimension>::type
    make(const Element &e) const {

      std::size_t size = referenceElement(e.geometry()).size(Element::mydimension);

      std::vector<typename IdSet::IdType> ids(size);
      for(std::size_t i = 0; i < size; ++i)
        ids[i] = idset.subId(e, i, Element::mydimension);
      return GeneralVertexOrder<Element::mydimension, Index>
               (e.type(), ids.begin(), ids.end());
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_VERTEXORDERFACTORY_HH
