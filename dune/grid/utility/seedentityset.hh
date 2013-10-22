// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_SEEDENTITYSET_HH
#define DUNE_GRID_UTILITY_SEEDENTITYSET_HH

#include <cassert>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/static_assert.hh>

namespace Dune {

  template<class Grid, int codim, class SeedIterator,
           class Category =
             typename std::iterator_traits<SeedIterator>::iterator_category>
  class SeedEntityIterator :
    public ForwardIteratorFacade<
    SeedEntityIterator<Grid, codim, SeedIterator, Category>,
      const typename Grid::template Codim<codim>::Entity,
      const typename Grid::template Codim<codim>::Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type>
  {
    typedef const typename Grid::template Codim<codim>::Entity Entity;
    typedef ForwardIteratorFacade<
      SeedEntityIterator, const Entity, const Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type
      > Base;
    friend Base;

    typedef typename Grid::template Codim<codim>::EntityPointer EntityPointer;
    const Grid *gridp_;
    SeedIterator seedIt_;
    mutable EntityPointer ep_;
    mutable bool valid_;

  public:
    typedef typename Base::Reference Reference;

    SeedEntityIterator(const Grid &grid, const SeedIterator &seedIt) :
      gridp_(&grid), seedIt_(seedIt), ep_(), valid_(false)
    { }

  private:
    Reference dereference() const
    {
      if(!valid_) {
        ep_ = gridp_->entityPointer(*seedIt_);
        valid_ = true;
      }
      return *ep_;
    }

    bool equals(const SeedEntityIterator &other) const
    {
      return seedIt_ == other.seedIt_;
    }

    void increment()
    {
      ++seedIt_;
      valid_ = false;
    }

  };

  template<class Grid, int codim, class SeedIterator>
  class SeedEntityIterator<Grid, codim, SeedIterator,
                           std::random_access_iterator_tag> :
    public RandomAccessIteratorFacade<
      SeedEntityIterator<Grid, codim, SeedIterator,
                         std::random_access_iterator_tag>,
      const typename Grid::template Codim<codim>::Entity,
      const typename Grid::template Codim<codim>::Entity &,
      typename std::iterator_traits<SeedIterator>::difference_type>
  {
  public:
    typedef typename Grid::template Codim<codim>::Entity Entity;

  private:
    typedef RandomAccessIteratorFacade<
      SeedEntityIterator, const Entity, const Entity &,
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

    SeedEntityIterator(const Grid &grid, const SeedIterator &seedIt) :
      gridp_(&grid), seedIt_(seedIt),
      ep_(gridp_->levelView(0).template end<codim>()), valid_(false)
    { }

    Reference dereference() const
    {
      if(!valid_) {
        ep_ = gridp_->entityPointer(*seedIt_);
        valid_ = true;
      }
      return *ep_;
    }

    Reference elementAt(DifferenceType n) const
    {
      dune_static_assert(AlwaysFalse<Grid>::value, "Subscription not "
                         "implemented for random-access SeedEntityIterator");
    }

    bool equals(const SeedEntityIterator &other) const
    {
      return seedIt_ == other.seedIt_;
    }

    void increment()
    {
      ++seedIt_;
      valid_ = false;
    }

    void decrement()
    {
      --seedIt_;
      valid_ = false;
    }

    void advance(DifferenceType n)
    {
      std::advance(seedIt_, n);
      valid_ = false;
    }

    DifferenceType distanceTo(const SeedEntityIterator &other) const
    {
      assert(grid_ == other.grid_);
      return std::distance(seedIt_, other.seedIt_);
    }
  };

  template<class EntityIterator,
           class Category =
             typename std::iterator_traits<EntityIterator>::iterator_category>
  class EntitySeedIteratorAdapter :
    public ForwardIteratorFacade<
      EntitySeedIteratorAdapter<EntityIterator, Category>,
      const typename EntityIterator::Entity::EntitySeed,
      typename EntityIterator::Entity::EntitySeed,
      typename std::iterator_traits<EntityIterator>::difference_type>
  {
    typedef typename EntityIterator::Entity::EntitySeed Seed;
    typedef ForwardIteratorFacade<
      EntitySeedIteratorAdapter, const Seed, Seed,
      typename std::iterator_traits<EntityIterator>::difference_type
      > Base;

    EntityIterator eIt_;
  public:
    typedef typename Base::Reference Reference;

    EntitySeedIteratorAdapter(const EntityIterator &eIt) :
      eIt_(eIt)
    { }

    Seed dereference() const
    {
      return eIt_->seed();
    }

    bool equals(const EntitySeedIteratorAdapter &other) const
    {
      return eIt_ == other.eIt_;
    }

    void increment()
    {
      ++eIt_;
    }

  };

  template<class EntityIterator>
  EntitySeedIteratorAdapter<EntityIterator>
  makeSeedIteratorAdaptor(const EntityIterator &eIt)
  {
    return EntitySeedIteratorAdapter<EntityIterator>(eIt);
  }

  template<class I>
  class IteratorRangeEntitySet {
    I begin_;
    I end_;

  public:
    typedef I const_iterator;
    typedef typename I::Entity Entity;

    IteratorRangeEntitySet(const I &begin, const I &end) :
      begin_(begin), end_(end)
    { }

    const I &begin() const
    {
      return begin_;
    }
    const I &end() const
    {
      return end_;
    }
  };

  template<class Grid, int codim>
  class SeedListPartitioning
  {
    typedef typename Grid::template Codim<codim>::EntitySeed Seed;
    typedef SeedEntityIterator<Grid, codim,
                               typename std::vector<Seed>::const_iterator>
      Iterator;

    const Grid *gridp_;
    std::vector<Seed> seedLists_;
    std::vector<std::size_t> pBegin_;
    std::vector<std::size_t> pEnd_;

  public:
    typedef IteratorRangeEntitySet<Iterator> EntitySet;

    template<class GV>
    SeedListPartitioning(const GV &gv) :
      gridp_(&gv.grid()),
      seedLists_(makeSeedIteratorAdaptor(gv.template begin<codim>()),
                 makeSeedIteratorAdaptor(gv.template end<codim>())),
      pBegin_(1, 0), pEnd_(1, seedLists_.size())
    { }

    std::size_t partitions() const
    {
      return pBegin_.size();
    }
    // return number of elements visited by partition
    std::size_t size(std::size_t partition) const
    {
      return pEnd_[partition] - pBegin_[partition];
    }
    // entityset for a particular partition
    EntitySet entitySet(std::size_t partition) const
    {
      return EntitySet(Iterator(*gridp_,
                                seedLists_.begin() + pBegin_[partition]),
                       Iterator(*gridp_,
                                seedLists_.begin() + pEnd_[partition]));
    }

    //! split a partition
    /**
     * \param partition    Numer of partition to replace.
     * \param newSeedLists Container of containers with seeds.  One
     *                     subcontainer for each thread.  The total number of
     *                     seeds in all subcontainer better be equal to the
     *                     size of the old partition.
     *
     * \returns Vector of new partition numbers.  The first entry will be
     *          equal to \c partition, and the number of entries will be equal
     *          to \c seedLists.size().
     */
    template<class SeedLists>
    std::vector<std::size_t> splitPartition(std::size_t partition,
                                            const SeedLists &newSeedLists)
    {
      std::size_t npartitions = std::distance(newSeedLists.begin(),
                                              newSeedLists.end());
      std::vector<std::size_t> newPartitions(npartitions);
      newPartitions[0] = partition;
      for(std::size_t i = 1; i < npartitions; ++i)
        newPartitions[i] = pBegin_.size() + i - 1;
      pBegin_.resize(pBegin_.size() + npartitions - 1);
      pEnd_.resize(pEnd_.size() + npartitions - 1);

      std::size_t pos = pBegin_[partition];
      std::size_t p = 0;
      for(const auto &seedList : newSeedLists)
      {
        pBegin_[newPartitions[p]] = pos;
        for(auto seed : seedList)
          seedLists_[pos++] = seed;
        pEnd_[newPartitions[p]] = pos;
        ++p;
      }

      return std::move(newPartitions);
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_SEEDENTITYSET_HH
