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

    //! construct
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
      for(const auto &e : entities<codim>(gv))
        ++pEnd_[partitioner.partition(e)];
      std::partial_sum(pEnd_.begin(), pEnd_.end(), pEnd_.begin());

      // construct begin indices
      pBegin_[0] = 0;
      std::copy(pEnd_.begin(), pEnd_.end()-1, pBegin_.begin()+1);

      // Vector of running indices, one for each partition
      auto pIndex = pBegin_;
      // fill the seed lists with the correct seeds
      for(const auto &e : entities<codim>(gv))
      {
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
     * \param pId          Numer of partition to replace.
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
    std::vector<Size> splitPartition(Size pId, const SeedLists &newSeedLists)
    {
      Size npartitions = std::distance(newSeedLists.begin(),
                                       newSeedLists.end());
      std::vector<Size> newPartitions(npartitions);
      newPartitions[0] = pId;
      for(Size i = 1; i < npartitions; ++i)
        newPartitions[i] = pBegin_.size() + i - 1;
      pBegin_.resize(pBegin_.size() + npartitions - 1);
      pEnd_.resize(pEnd_.size() + npartitions - 1);

      Size pos = pBegin_[pId];
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

  private:
    const Grid *gridp_;
    std::vector<Seed> seedLists_;
    std::vector<Size> pBegin_;
    std::vector<Size> pEnd_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_SEEDLIST_HH
