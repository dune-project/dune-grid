// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ITERATORADAPTERS_HH
#define DUNE_GRID_UTILITY_ITERATORADAPTERS_HH

#include <cassert>
#include <iterator>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/typetraits.hh>

namespace Dune {

  //! generic filtering iterator for entity sets
  template<class EntitySet, class HostIterator>
  class EntitySetIterator :
    public ForwardIteratorFacade<
      EntitySetIterator<EntitySet, HostIterator>,
      const typename std::iterator_traits<HostIterator>::value_type,
      const typename std::iterator_traits<HostIterator>::value_type &,
      typename std::iterator_traits<HostIterator>::difference_type>
  {
  public:
    //! type of entity
    typedef typename remove_const<
      typename std::iterator_traits<HostIterator>::value_type>::type Entity;
    //! codimension of entity
    enum { codimension = Entity::codimension };

  private:
    typedef ForwardIteratorFacade<
      EntitySetIterator, const Entity, const Entity &,
      typename std::iterator_traits<HostIterator>::difference_type> Base;

  public:
    typedef typename Base::Reference Reference;

    //! construct
    EntitySetIterator(const EntitySet &entitySet, const HostIterator &host,
                      const HostIterator &hostEnd) :
      entitySet_(entitySet), host_(host), hostEnd_(hostEnd)
    {
      findValid();
    }

    //! construct end-iterator
    EntitySetIterator(const EntitySet &entitySet,
                      const HostIterator &hostEnd) :
      entitySet_(entitySet), host_(hostEnd), hostEnd_(hostEnd)
    { }

    Reference dereference() const
    {
      return *host_;
    }

    bool equals(const EntitySetIterator &rhs) const
    {
      return host_ == rhs.host_;
    }

    void increment() {
      assert(host_ != hostEnd_);
      ++host_;
      findValid();
    };

  private:
    void findValid()
    {
      while(host_ != hostEnd_ && !entitySet_.contains(*host_))
        ++host_;
    }

    EntitySet entitySet_;
    HostIterator host_;
    HostIterator hostEnd_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ITERATORADAPTERS_HH
