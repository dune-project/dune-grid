// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_SEEDENTITYSET_HH
#define DUNE_GRID_UTILITY_SEEDENTITYSET_HH

#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include <dune/grid/utility/entityrange.hh>
#include <dune/grid/utility/iteratoradapters.hh>

namespace Dune {

  template<class Grid, int codim>
  class SeedListPartitioning
  {
    typedef typename Grid::template Codim<codim>::EntitySeed Seed;
    typedef SeedToEntityIteratorAdapter<
      Grid, typename std::vector<Seed>::const_iterator> Iterator;

    const Grid *gridp_;
    std::vector<Seed> seedLists_;
    std::vector<std::size_t> pBegin_;
    std::vector<std::size_t> pEnd_;

  public:
    typedef IteratorEntityRange<Iterator> EntitySet;

    template<class GV>
    SeedListPartitioning(const GV &gv) :
      gridp_(&gv.grid()),
      seedLists_(makeEntityToSeedIteratorAdaptor(gv.template begin<codim>()),
                 makeEntityToSeedIteratorAdaptor(gv.template end<codim>())),
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
