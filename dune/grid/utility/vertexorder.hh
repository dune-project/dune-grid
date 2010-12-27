// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_VERTEXORDER_HH
#define DUNE_GRID_UTILITY_VERTEXORDER_HH

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <vector>

#include <dune/common/geometrytype.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune {

  //! algorithm to reduce vertex order information
  /**
   * \code
   *#include <dune/grid/utility/vertexorder.hh>
   * \endcode
   *
   * \param inBegin Start of the range of ids to reduce.
   * \param inEnd   End of the range of ids to reduce.
   * \param outIt   Start of the sequence where to store the result.
   *
   * \c inBegin and \c inEnd must be ForwardIterators; their \c value_type may
   * constant.  \c outIt must be an OutputIterator and must allow \c
   * std::distance(inBegin,inEnd) increments.  Complexity is quadratic.
   *
   * \sa GeneralVertexOrder, VertexOrderByIdFactory
   */
  template<class InIterator, class OutIterator>
  void reduceOrder(const InIterator& inBegin, const InIterator& inEnd,
                   OutIterator outIt)
  {
    typedef std::less<
        typename std::iterator_traits<InIterator>::value_type
        > less_t;
    static const less_t less = less_t();

    for(InIterator inIt = inBegin; inIt != inEnd; ++inIt, ++outIt)
      *outIt = std::count_if(inBegin, inEnd, std::bind2nd(less, *inIt));
  }

  //! Class providing information on the ordering of vertices
  /**
   * \tparam dim    Dimension of the entity this class provides ordering
   *                information for.
   * \tparam Index_ Type of the indices.  Must be integral, may be
   *                non-negative.
   *
   * This class provides ordering information for all codimensions, including
   * the element itself.
   *
   * \warning The Interface of the VertexOrder stuff is subject to change.  It
   *          is currently needed to use some global-valued finite elements
   *          from dune-localfunctions.
   *
   * \sa reduceOrder(), VertexOrderByIdFactory
   */
  template<std::size_t dim, class Index_ = std::size_t>
  class GeneralVertexOrder {
    typedef GenericReferenceElement<double, dim> RefElem;
    typedef GenericReferenceElements<double, dim> RefElems;

    const RefElem& refelem;
    GeometryType gt;
    std::vector<Index_> vertexOrder;

  public:
    //! Type of indices
    typedef Index_ Index;

    //! Iterate over the vertex indices of some sub-entity
    class iterator;

    //! export the dimension of the entity we provide information for
    static const std::size_t dimension = dim;
    //! get type of the entity's geometry
    const GeometryType &type() const { return gt; }

    //! construct a GeneralVertexOrder
    /**
     * \param gt_     Geometry type of the entity we provide information for.
     * \param inBegin Start of the range of vertex ids.
     * \param inEnd   End of the range of vertex ids.
     *
     * \c inBegin and \c inEnd denote the range of vertes ids to provide.
     * This class stores a reduced copy of the ids, converted to type Index.
     */
    template<class InIterator>
    GeneralVertexOrder(const GeometryType& gt_, const InIterator &inBegin,
                       const InIterator &inEnd) :
      refelem(RefElems::general(gt_)), gt(gt_),
      vertexOrder(refelem.size(dim))
    { reduceOrder(inBegin, inEnd, vertexOrder.begin()); }

    //! get begin iterator for the vertex indices of some sub-entity
    /**
     * \param codim     Codimension of the sub-entity.
     * \param subEntity Index of the sub-entity within that codimension.
     */
    iterator begin(std::size_t codim, std::size_t subEntity) const
    { return iterator(*this, codim, subEntity); }
    //! get end iterator for the vertex indices of some sub-entity
    /**
     * \param codim     Codimension of the sub-entity.
     * \param subEntity Index of the sub-entity within that codimension.
     */
    iterator end(std::size_t codim, std::size_t subEntity) const {
      return iterator(*this, codim, subEntity,
                      refelem.size(subEntity, codim, dim));
    }

    //! get a vector of reduced indices for some sub-entity
    /**
     * \param codim     Codimension of the sub-entity.
     * \param subEntity Index of the sub-entity within that codimension.

     * \param order     Where to store the result.  This function resizes the
     *                  vector to the suitable size.
     */
    void getReduced(std::size_t codim, std::size_t subEntity,
                    std::vector<Index>& order) const
    {
      order.resize(refelem.size(subEntity, codim, dim));
      reduceOrder(begin(codim, subEntity), end(codim, subEntity),
                  order.begin());
    }
  };

  //! Iterate over the vertex indices of some sub-entity
  /**
   * This is a random access iterator with constant \c value_type.
   */
  template<std::size_t dim, class Index_>
  class GeneralVertexOrder<dim, Index_>::iterator :
    public Dune::RandomAccessIteratorFacade<iterator, const Index_>
  {
    const GeneralVertexOrder *order;
    std::size_t codim;
    std::size_t subEntity;
    std::size_t vertex;

    iterator(const GeneralVertexOrder &order_, std::size_t codim_,
             std::size_t subEntity_, std::size_t vertex_ = 0) :
      order(&order_), codim(codim_), subEntity(subEntity_), vertex(vertex_)
    { }

  public:
    const Index &dereference() const {
      return order->vertexOrder[order->refelem.subEntity(subEntity, codim,
                                                         vertex, dim)];
    }
    const Index &elementAt(std::ptrdiff_t n) const {
      return order->vertexOrder[order->refelem.subEntity(subEntity, codim,
                                                         vertex+n, dim)];
    }
    bool equals(const iterator &other) const {
      return order == other.order && codim == other.codim &&
             subEntity == other.subEntity && vertex == other.vertex;
    }
    void increment() { ++vertex; }
    void decrement() { --vertex; }
    void advance(std::ptrdiff_t n) { vertex += n; }
    std::ptrdiff_t distanceTo(const iterator &other) const {
      // make sure we reference the same container
      assert(order == other.order && codim == other.codim &&
             subEntity == other.subEntity);
      if(vertex < other.vertex) return other.vertex - vertex;
      else return -static_cast<std::ptrdiff_t>(vertex - other.vertex);
    }

    friend class GeneralVertexOrder<dim, Index>;

    //! public default constructor
    /**
     * The contructed iterator object will have a singular value.  The only
     * valid operations will be assignment of a non-singular value and
     * destruction, all other operations will result in undefined behaviour.
     */
    iterator() { }
  };

  //! Factory for GeneralVertexOrder objects using an IdSet
  /**
   * \tparam IdSet Type IdSet used to get the ids of the vertices.
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
     * becomes invalid.  The only valid operation on an factory with singular
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
      typedef GenericReferenceElements<
          typename Element::ctype,
          Element::mydimension
          > RefElems;
      std::size_t size =
        RefElems::general(e.type()).size(Element::mydimension);

      std::vector<typename IdSet::IdType> ids(size);
      for(std::size_t i = 0; i < size; ++i)
        ids[i] = idset.subId(e, i, Element::mydimension);
      return GeneralVertexOrder<Element::mydimension, Index>
               (e.type(), ids.begin(), ids.end());
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_VERTEXORDER_HH
