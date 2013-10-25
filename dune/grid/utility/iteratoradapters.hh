// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ITERATORADAPTERS_HH
#define DUNE_GRID_UTILITY_ITERATORADAPTERS_HH

#include <cassert>
#include <iterator>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune {

  //////////////////////////////////////////////////////////////////////
  //
  //  EntitySetIterator
  //

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

  //////////////////////////////////////////////////////////////////////
  //
  //  SeedToEntityIteratorAdapter
  //

  //! adapt an iterator over entityseeds into an interator over entities
  /**
   * \tparam Grid         Type of grid.
   * \tparam SeedIterator Type of the iterator providing the seeds.
   * \tparam Category     Category of seed iterator.  This is used to
   *                      specialize for random access iterators.  It is best
   *                      left at the default.
   */
  template<class Grid, class SeedIterator,
           class Category =
             typename std::iterator_traits<SeedIterator>::iterator_category>
  class SeedToEntityIteratorAdapter :
    public ForwardIteratorFacade<
      SeedToEntityIteratorAdapter<Grid, SeedIterator, Category>,
      const typename Grid::template Codim<
        std::iterator_traits<SeedIterator>::value_type::codimension>::Entity,
      const typename Grid::template Codim<
        std::iterator_traits<SeedIterator>::value_type::codimension>::Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type>
  {
    typedef typename std::iterator_traits<SeedIterator>::value_type Seed;
    enum { codim = Seed::codimension };
    typedef const typename Grid::template Codim<codim>::Entity Entity;
    typedef ForwardIteratorFacade<
      SeedToEntityIteratorAdapter, const Entity, const Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type
      > Base;

    typedef typename Grid::template Codim<codim>::EntityPointer EntityPointer;
    const Grid *gridp_;
    SeedIterator seedIt_;
    mutable EntityPointer ep_;
    mutable bool valid_;

  public:
    typedef typename Base::Reference Reference;

    //! construct
    SeedToEntityIteratorAdapter(const Grid &grid, const SeedIterator &seedIt) :
      gridp_(&grid), seedIt_(seedIt), ep_(), valid_(false)
    { }

    //! For Iteratorfacade
    Reference dereference() const
    {
      if(!valid_) {
        ep_ = gridp_->entityPointer(*seedIt_);
        valid_ = true;
      }
      return *ep_;
    }

    //! For Iteratorfacade
    bool equals(const SeedToEntityIteratorAdapter &other) const
    {
      return seedIt_ == other.seedIt_;
    }

    void increment()
    {
      ++seedIt_;
      valid_ = false;
    }
  };

  //! adapt an iterator over entityseeds into an interator over entities
  /**
   * This is the specialization for random-access iterators.
   */
  template<class Grid, class SeedIterator>
  class SeedToEntityIteratorAdapter<Grid, SeedIterator,
                                    std::random_access_iterator_tag> :
    public RandomAccessIteratorFacade<
      SeedToEntityIteratorAdapter<Grid, SeedIterator,
                                  std::random_access_iterator_tag>,
      const typename Grid::template Codim<
        std::iterator_traits<SeedIterator>::value_type::codimension>::Entity,
      const typename Grid::template Codim<
        std::iterator_traits<SeedIterator>::value_type::codimension>::Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type>
  {
    typedef typename std::iterator_traits<SeedIterator>::value_type Seed;
    enum { codim = Seed::codimension };
    typedef typename Grid::template Codim<codim>::Entity Entity;
    typedef RandomAccessIteratorFacade<
      SeedToEntityIteratorAdapter, const Entity, const Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type
      > Base;

    typedef typename Grid::template Codim<codim>::EntityPointer EntityPointer;
    const Grid *gridp_;
    SeedIterator seedIt_;
    mutable EntityPointer ep_;
    mutable bool valid_;

  public:
    typedef typename Base::Reference Reference;
    typedef typename Base::DifferenceType DifferenceType;

    //! construct
    SeedToEntityIteratorAdapter(const Grid &grid, const SeedIterator &seedIt) :
      gridp_(&grid), seedIt_(seedIt),
      ep_(gridp_->levelView(0).template end<codim>()), valid_(false)
    { }

    //! For Iteratorfacade
    Reference dereference() const
    {
      if(!valid_) {
        ep_ = gridp_->entityPointer(*seedIt_);
        valid_ = true;
      }
      return *ep_;
    }

    //! For Iteratorfacade
    Reference elementAt(DifferenceType n) const
    {
      dune_static_assert(AlwaysFalse<Grid>::value, "Subscription not "
                         "implemented for random-access SeedEntityIterator");
    }

    //! For Iteratorfacade
    bool equals(const SeedToEntityIteratorAdapter &other) const
    {
      return seedIt_ == other.seedIt_;
    }

    //! For Iteratorfacade
    void increment()
    {
      ++seedIt_;
      valid_ = false;
    }

    //! For Iteratorfacade
    void decrement()
    {
      --seedIt_;
      valid_ = false;
    }

    //! For Iteratorfacade
    void advance(DifferenceType n)
    {
      std::advance(seedIt_, n);
      valid_ = false;
    }

    //! For Iteratorfacade
    DifferenceType distanceTo(const SeedToEntityIteratorAdapter &other) const
    {
      assert(grid_ == other.grid_);
      return std::distance(seedIt_, other.seedIt_);
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ITERATORADAPTERS_HH
