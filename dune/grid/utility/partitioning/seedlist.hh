// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONING_SEEDLIST_HH
#define DUNE_GRID_UTILITY_PARTITIONING_SEEDLIST_HH

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include <dune/grid/utility/entityrange.hh>
#include <dune/grid/utility/iteratoradapters.hh>

namespace Dune {

  //! partitioning data structure based on lists of seeds for each partition
  template<class Grid, int codim>
  class SeedListPartitioning
  {
    typedef typename Grid::template Codim<codim>::EntitySeed Seed;
    typedef SeedToEntityIteratorAdapter<
      Grid, typename std::vector<Seed>::const_iterator> Iterator;

  public:
    //! type of partitions
    typedef IteratorEntityRange<Iterator> Partition;
    //! type used to count partitions
    typedef typename std::iterator_traits<
      typename std::vector<Seed>::iterator
      >::difference_type Size;

    //! construct
    template<class GV>
    SeedListPartitioning(const GV &gv) :
      gridp_(&gv.grid()),
      seedLists_(makeEntityToSeedIteratorAdaptor(gv.template begin<codim>()),
                 makeEntityToSeedIteratorAdaptor(gv.template end<codim>())),
      pBegin_(1, 0), pEnd_(1, seedLists_.size())
    { }

    //! construct from a partitioner
    /**
     * The partitioner should support the two functions \c partition(entity)
     * and \c partitions().  The reference to the partitioner is not stored
     * internally, so may be a temporary.
     */
    template<class GV, class Partitioner>
    SeedListPartitioning(const GV &gv, const Partitioner &partitioner) :
      gridp_(&gv.grid()),
      pBegin_(partitioner.partitions()), pEnd_(partitioner.partitions(), 0)
    {
      // Initialize seed lists.  We arbitrarily use the seed of the first
      // entity for this.  Make certain not to dereference gv.begin() when the
      // gridview is empty
      if(gv.size(codim))
        seedLists_.assign(gv.size(codim), gv.template begin<0>()->seed());

      // compute end indices
      auto erange = entityRange<codim>(gv);
      for(auto it = erange.begin(); it != erange.end(); ++it)
        ++pEnd_[partitioner.partition(*it)];
      std::partial_sum(pEnd_.begin(), pEnd_.end(), pEnd_.begin());

      // construct begin indices
      pBegin_[0] = 0;
      std::copy(pEnd_.begin(), pEnd_.end()-1, pBegin_.begin()+1);

      // Vector of running indices, one for each partition
      auto pIndex = pBegin_;
      // fill the seed lists with the correct seeds
      for(auto it = erange.begin(); it != erange.end(); ++it)
      {
        const auto &e = *it;
        auto &index = pIndex[partitioner.partition(e)];
        seedLists_[index] = e.seed();
        ++index;
      }
    }

    //! return maximum number of partitions
    Size partitions() const
    {
      return pBegin_.size();
    }

    //! return a particular partition
    Partition partition(Size pId) const
    {
      return Partition(Iterator(*gridp_,
                                seedLists_.begin() + pBegin_[pId]),
                       Iterator(*gridp_,
                                seedLists_.begin() + pEnd_[pId]));
    }

    //! split a partition
    /**
     * \param pId         Numer of partition to replace.
     * \param partitioner A partitioner object that supports the methods \c
     *                    partition(entity) and \c partitions().
     *
     * \returns Vector of new partition numbers.  The first entry will be
     *          equal to \c pId, and the number of entries will be equal to \c
     *          partitioner.partitions().
     */
    template<class Partitioner>
    std::vector<Size> splitPartition(Size pId, const Partitioner &partitioner)
    {
      // copy of the old partition
      std::vector<Seed> tmpLists(seedLists_.begin()+pBegin_[pId],
                                 seedLists_.begin()+pEnd_[pId]);

      // compute new partition numbers
      std::vector<Size> newPartitions(partitioner.partitions());
      newPartitions[0] = pId;
      for(std::size_t i = 1; i < newPartitions.size(); ++i)
        newPartitions[i] = pBegin_.size() + i - 1;

      // compute partition size
      std::vector<std::size_t> pSize(partitioner.partitions(), 0);
      for(const auto &seed : tmpLists)
        ++pSize[partitioner.partition(*gridp_->entityPointer(seed))];

      // set begin indices of new partition
      pBegin_.resize(pBegin_.size() + newPartitions.size() - 1);
      for(std::size_t i = 1; i < newPartitions.size(); ++i)
        pBegin_[newPartitions[i]] =
          pBegin_[newPartitions[i-1]] + pSize[i-1];
      // set end indices of new partition to begin indices...
      pEnd_.resize(pEnd_.size() + newPartitions.size() - 1);
      for(auto p : newPartitions)
        pEnd_[p] = pBegin_[p];

      // ...and use them as per-partition running variables when copying the
      // seed back
      for(const auto &seed : tmpLists) {
        auto &index = pEnd_[newPartitions[partitioner.partition
                                          (*gridp_->entityPointer(seed))]];
        seedLists_[index] = seed;
        ++index;
      }
      // now the end indices should be correct

      return std::move(newPartitions);
    }

  private:
    const Grid *gridp_;
    std::vector<Seed> seedLists_;
    std::vector<Size> pBegin_;
    std::vector<Size> pEnd_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_SEEDLIST_HH
